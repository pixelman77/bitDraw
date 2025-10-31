#include <SDL3/SDL.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <bits/stdc++.h>
#include <OBJ_Loader.h>

const float FLOAT_MAX = 100000.0f; //not the max, but a usable amount
const float ZDEPTH_VISIBLE = 0.4f;
const float ZDEPTH_BIAS = 0.0f;
const int orthographicHurlUnit = 1000; //camera position is moved back by this amount to simulate a far away camera 

enum ProjectionMode
{
    PROJ_PERSPECTIVE,
    PROJ_ORTHOGRAPHIC
};

enum RenderMode
{
    RENDER_WIREFRAME,
    RENDER_PAINTERS,
    RENDER_ZBUFFER_MODEL,
    RENDER_ZBUFFER_BUFFER,
    RENDER_NO_DEPTH
};

enum TriangleRasterMethod
{
    SCANLINE,
    BOUNDING_BOX
};

enum LightMode
{
    CONSTANT,
    PHONG,
    SIMPLE_LAMBERT
};

enum ColorMode
{
    SOLID_WHITE,
    RANDOM_COLOR,
    TEXTURE
};

enum ShadeMode
{
    FLAT,
    SMOOTH
};

//default options
ProjectionMode projectionMode = ProjectionMode::PROJ_PERSPECTIVE;
RenderMode renderMode = RenderMode::RENDER_ZBUFFER_MODEL;
TriangleRasterMethod rasterMode = TriangleRasterMethod::BOUNDING_BOX;
LightMode lightMode = LightMode::SIMPLE_LAMBERT;
ColorMode colorMode = ColorMode::SOLID_WHITE;
ShadeMode shadeMode = ShadeMode::SMOOTH;
bool backfaceCulling = false;



//buffer input gets a seizure with these and I have no idea why
uint32_t SCREEN_WIDTH = 800;
uint32_t SCREEN_HEIGHT = 600;
float HALFWIDTH = SCREEN_WIDTH / 2;
float HALFHEIGHT = SCREEN_HEIGHT / 2;

bool assignedRandomColor = false;

uint32_t packARGB(uint8_t a, uint8_t r, uint8_t g, uint8_t b) {
    return (static_cast<uint32_t>(a) << 24) |
           (static_cast<uint32_t>(r) << 16) |
           (static_cast<uint32_t>(g) << 8)  |
           (static_cast<uint32_t>(b));
}

double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}

double signed_triangle_area(int ax, int ay, int bx, int by, int cx, int cy) {
    return .5*((by-ay)*(bx+ax) + (cy-by)*(cx+bx) + (ay-cy)*(ax+cx));
}


struct ScreenBuffer
{
    // note to self: w of h/x of y
    // down growing
    std::vector<uint32_t> buffer;
    std::vector<double> zBuffer;
    int width, height;
    float halfwidth, halfheight;

    ScreenBuffer(uint32_t iwidth, uint32_t iheigth){

        this->width = iwidth;
        this->height = iheigth;
        halfheight = height / 2;
        halfwidth = width / 2;
        buffer = std::vector<uint32_t>(width * height, 0);
        zBuffer = std::vector<double>(width * height, FLOAT_MAX);
    }

    void fill(uint32_t color){
        std::fill(buffer.begin(), buffer.end(), color);
    }

    //boundry check too slow?
    uint32_t get(int x, int y){
        if(x >= width || y>=height) return 0;
        return buffer[y * width + x];
    }

    void set(int x, int y, uint32_t color){
        if(x < 0 || y < 0 || x >= width || y>=height) return;
        buffer[y * width + x] = color;
    }

    void drainZBuffer(){
        std::fill(zBuffer.begin(), zBuffer.end(), FLOAT_MAX);
    }

    //boundry check too slow?
    double getZBuffer(int x, int y){
        if(x >= width || y>=height) return 0;
        return zBuffer[y * width + x];
    }

    void setZBuffer(int x, int y, float val){
        if(x < 0 || y < 0 || x >= width || y>=height) return;
        zBuffer[y * width + x] = val;
    }

    bool updateZBuffer(int x, int y, double val){
        if(getZBuffer(x, y) <= val){ return false; }
        setZBuffer(x, y, val);
        return true;
    }

//this following code only draws vertexes (either from one side or wherever), will I ever investigate?
/*
    void line(int x1, int y1, int x2, int y2, uint32_t color){
        //shamelessly stolen from tinyRenderer
        bool steep = std::abs(x1-x2) < std::abs(y1-y2);
        if (steep) { // if the line is steep, we transpose the image
            std::swap(x1, y1);
            std::swap(x2, y2);
        }
        if (x1>x2) { // make it left−to−right
            std::swap(x1, x2);
            std::swap(y1, y2);
        }
        int y = y1;
        int ierror = 0;
        for (int x=x1; x<=x1; x++) {
            if (steep) // if transposed, de−transpose
                set(y, x, color);
            else
                set(x, y, color);
            ierror += 2 * std::abs(y2-y1);
            y += (y2 > y1 ? 1 : -1) * (ierror > x2 - x1);
            ierror -= 2 * (x2-x1)   * (ierror > x2 - x1);
        } 

    }
    */
};


void drawLine( int x0, int y0, int x1, int y1, ScreenBuffer &buffer, uint32_t color) {
        int dx = std::abs(x1 - x0);
        int dy = std::abs(y1 - y0);
        int sx = (x0 < x1) ? 1 : -1;
        int sy = (y0 < y1) ? 1 : -1;
        int err = dx - dy;

        while (true) {
            buffer.set(x0, y0, color); // draw the pixel

            if (x0 == x1 && y0 == y1) break;

            int e2 = 2 * err;
            if (e2 > -dy) {
                err -= dy;
                x0 += sx;
            }
            if (e2 < dx) {
                err += dx;
                y0 += sy;
            }
        }
}

//currently takes no color from shader
//do fix if you plan on showcasing
void triangleRasterScanLine(int ax, int ay, int bx, int by, int cx, int cy, ScreenBuffer &buffer, 
    uint32_t (*color)(float, float, float) ) {
        // sort the vertices, a,b,c in ascending y order (bubblesort yay!)
        if (ay>by) { std::swap(ax, bx); std::swap(ay, by); }
        if (ay>cy) { std::swap(ax, cx); std::swap(ay, cy); }
        if (by>cy) { std::swap(bx, cx); std::swap(by, cy); }
        int total_height = cy-ay;

        if (ay != by) { // if the bottom half is not degenerate
            int segment_height = by - ay;
            for (int y=ay; y<=by; y++) { // sweep the horizontal line from ay to by
                int x1 = ax + ((cx - ax)*(y - ay)) / total_height;
                int x2 = ax + ((bx - ax)*(y - ay)) / segment_height;
                for (int x=std::min(x1,x2); x<std::max(x1,x2); x++)  // draw a horizontal line
                    buffer.set(x, y, 0);
            }
        }
        if (by != cy) { // if the upper half is not degenerate
            int segment_height = cy - by;
            for (int y=by; y<=cy; y++) { // sweep the horizontal line from by to cy
                int x1 = ax + ((cx - ax)*(y - ay)) / total_height;
                int x2 = bx + ((cx - bx)*(y - by)) / segment_height;
                for (int x=std::min(x1,x2); x<std::max(x1,x2); x++)  // draw a horizontal line
                    buffer.set(x, y, 0);
            }
        }
}

//not working
void triangleRasterScanline_ZBuffered(int ax, int ay, int bx, int by, int cx, int cy, float az, float bz, float cz, ScreenBuffer &buffer, 
    uint32_t (*color)(float, float, float)) {
        
        double total_area = signed_triangle_area(ax, ay, bx, by, cx, cy);

        if (total_area<1) return;
        

        if (ay>by) { std::swap(ax, bx); std::swap(ay, by); std::swap(az, bz); }
        if (ay>cy) { std::swap(ax, cx); std::swap(ay, cy); std::swap(az, cz); }
        if (by>cy) { std::swap(bx, cx); std::swap(by, cy); std::swap(bz, cz); }
        int total_height = cy-ay;
        

        if (ay != by) { 
            int segment_height = by - ay;
            for (int y=ay; y<=by; y++) { 
                int x1 = ax + ((cx - ax)*(y - ay)) / total_height;
                int x2 = ax + ((bx - ax)*(y - ay)) / segment_height;
                for (int x=std::min(x1,x2); x<std::max(x1,x2); x++){  
                    double alpha = signed_triangle_area(x, y, bx, by, cx, cy) / total_area;
                    double beta  = signed_triangle_area(x, y, cx, cy, ax, ay) / total_area;
                    double gamma = signed_triangle_area(x, y, ax, ay, bx, by) / total_area;

                    double z = (alpha * az + beta * bz + gamma * cz);
                    if(buffer.updateZBuffer(x, y, z)){ buffer.set(x, y, color(alpha, beta, gamma)); }
                }                    
            }
        }
        if (by != cy) { 
            int segment_height = cy - by;
            for (int y=by; y<=cy; y++) { 
                int x1 = ax + ((cx - ax)*(y - ay)) / total_height;
                int x2 = bx + ((cx - bx)*(y - by)) / segment_height;
                for (int x=std::min(x1,x2); x<std::max(x1,x2); x++){  
                    double alpha = signed_triangle_area(x, y, bx, by, cx, cy) / total_area;
                    double beta  = signed_triangle_area(x, y, cx, cy, ax, ay) / total_area;
                    double gamma = signed_triangle_area(x, y, ax, ay, bx, by) / total_area;

                    double z = (alpha * az + beta * bz + gamma * cz);
                    if(buffer.updateZBuffer(x, y, z)){ buffer.set(x, y, color(alpha, beta, gamma)); } 
                }                   
            }
        }
    }



void triangleRasterBoundingBox(int ax, int ay, int bx, int by, int cx, int cy, ScreenBuffer &buffer, 
    uint32_t (*color)(float, float, float)) {
        int bbminx = std::min(std::min(ax, bx), cx); // bounding box for the triangle
        int bbminy = std::min(std::min(ay, by), cy); // defined by its top left and bottom right corners
        int bbmaxx = std::max(std::max(ax, bx), cx);
        int bbmaxy = std::max(std::max(ay, by), cy);
        double total_area = signed_triangle_area(ax, ay, bx, by, cx, cy);

    //he got tricked into thinking it would work this time!
    //mint issue?
        #pragma omp parallel for
        for (int x=bbminx; x<=bbmaxx; x++) {
            for (int y=bbminy; y<=bbmaxy; y++) {
                double alpha = signed_triangle_area(x, y, bx, by, cx, cy) / total_area;
                double beta  = signed_triangle_area(x, y, cx, cy, ax, ay) / total_area;
                double gamma = signed_triangle_area(x, y, ax, ay, bx, by) / total_area;
                if (alpha<0 || beta<0 || gamma<0) continue; // negative barycentric coordinate => the pixel is outside the triangle
                buffer.set(x, y, color(alpha, beta, gamma));
            }
        }
    }


void triangleRasterBoundingBox_ZBuffered(int ax, int ay, int bx, int by, int cx, int cy, float az, float bz, float cz, ScreenBuffer &buffer, 
    uint32_t (*color)(float, float, float)) {
        int bbminx = std::min(std::min(ax, bx), cx); // bounding box for the triangle
        int bbminy = std::min(std::min(ay, by), cy); // defined by its top left and bottom right corners
        int bbmaxx = std::max(std::max(ax, bx), cx);
        int bbmaxy = std::max(std::max(ay, by), cy);
        double total_area = signed_triangle_area(ax, ay, bx, by, cx, cy);

        if(total_area<1) return; //IMPORTANT!! tiny triangles cause issues with z-buffer calculations, making them appear as lines through bigger faces

        #pragma omp parallel for
        for (int x=bbminx; x<=bbmaxx; x++) {
            for (int y=bbminy; y<=bbmaxy; y++) {
                double alpha = signed_triangle_area(x, y, bx, by, cx, cy) / total_area;
                double beta  = signed_triangle_area(x, y, cx, cy, ax, ay) / total_area;
                double gamma = signed_triangle_area(x, y, ax, ay, bx, by) / total_area;
                if (alpha<0.f || beta<0.f || gamma<0.f) continue; // negative barycentric coordinate => the pixel is outside the triangle
                double z = 1/(alpha * az + beta * bz + gamma * cz);
                if(buffer.updateZBuffer(x, y, z)) buffer.set(x, y, color(alpha, beta, gamma));
            }
        }
    }


/*
struct Fragment
{
    ScreenPoint p1, p2, p3;
    std::vector<uint32_t> color;
    Fragment(ScreenPoint po1, ScreenPoint po2, ScreenPoint po3){
        p1 = po1;
        p2 = po2;
        p3 = po3;
    
        int bbminx = std::min(std::min(p1.x, p2.x), p3.x);
        int bbminy = std::min(std::min(p1.y, p2.y), p3.y);
        int bbmaxx = std::max(std::max(p1.x, p2.x), p3.x);
        int bbmaxy = std::max(std::max(p1.y, p2.y), p3.y);
    
    
    }
};
*/


struct Vector2d
{
    float x, y;
    Vector2d(float x = 0, float y = 0)
    {
        this->x = x;
        this->y = y;
    }
    // arithmetc logic?
};

struct ScreenPoint
{
    int x, y;
    double z;
    ScreenPoint(int x = 0, int y = 0, double z = 0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        
    }
};

struct Vector3d
{
    float x, y, z;

    Vector3d(float x = 0, float y = 0, float z = 0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // --- Operators ---
    Vector3d operator+(const Vector3d &rhs) const
    {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    Vector3d operator-(const Vector3d &rhs) const
    {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }

    Vector3d operator*(float rhs) const
    {
        return {x * rhs, y * rhs, z * rhs};
    }

    Vector3d operator/(float rhs) const
    {
        return {x / rhs, y / rhs, z / rhs};
    }

    // --- Vector math ---
    Vector3d normalized() const
    {
        float len = length();
        if (len == 0) return {0, 0, 0};
        return *this / len;
    }

    float squaredDistance(const Vector3d &rhs) const
    {
        float dx = x - rhs.x;
        float dy = y - rhs.y;
        float dz = z - rhs.z;
        return dx * dx + dy * dy + dz * dz;
    }

    Vector3d cross(const Vector3d &b) const
    {
        return {
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x};
    }

    float dot(const Vector3d &b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }

    float length() const
    {
        return std::sqrt(dot(*this));
    }

    Vector3d rotateByAngle(const Vector3d &axis, float angle) const
    {
        Vector3d u = axis.normalized();
        float cosA = std::cos(angle);
        float sinA = std::sin(angle);

        return (*this * cosA) +
               (u.cross(*this) * sinA) +
               (u * (u.dot(*this) * (1 - cosA)));
    }


};


struct Triangle
{
    Vector3d v1, v2, v3, vn1, vn2, vn3;
    uint8_t a, r, g, b;

    Triangle()
    {
        v1 = Vector3d();
        v2 = Vector3d();
        v3 = Vector3d();
        a = 255;
        r = SDL_rand(256);
        g = SDL_rand(256);
        b = SDL_rand(256);
        vn1 = Vector3d();
        vn2 = Vector3d();
        vn3 = Vector3d();
    }

    Vector3d centroid() {
        return Vector3d(
            (v1.x + v2.x + v3.x) / 3.0f,
            (v1.y + v2.y + v3.y) / 3.0f,
            (v1.z + v2.z + v3.z) / 3.0f);
    }

    Vector3d Normal() {
        Vector3d edge1 = v2 - v1;
        Vector3d edge2 = v3 - v1;
        return edge1.cross(edge2).normalized();
    }


};




struct Camera
{
    Vector3d position, rotation;
    float fov;    // for perspective
    float zoom;   // for zoom
    float aspect; // Aspect ratio (width / height)
    float near;
    float far;

    Camera()
    {
        position = Vector3d();
        rotation = Vector3d();
        fov = 1.0f;
        zoom = 80.0f;
        aspect = SCREEN_WIDTH / SCREEN_HEIGHT;
        near = 0.0f;
        far = 10.0f;
    }

    Vector3d getDirection() {
        float cosPitch = std::cos(rotation.x);
        float sinPitch = std::sin(rotation.x);
        float cosYaw   = std::cos(rotation.y);
        float sinYaw   = std::sin(rotation.y);

        return {
            cosPitch * sinYaw,  // x
            -sinPitch,          // y (note: depends on convention)
            cosPitch * cosYaw   // z
        };
        //is this normalized?
    }

    void orbitAroundPoint(const Vector3d& target, float angle, char axis, bool keepFacing) {
    //funnily enough, orbiting x causes a sin motion which while not intended, is pretty useful


    // Translate to target-relative coordinates
    Vector3d relative = position - target;

    // Rotate around chosen axis
    float c = std::cos(angle);
    float s = std::sin(angle);
    Vector3d rotated = relative;

    switch (axis) {
        case 'x': // rotate around X-axis
            rotated.y = relative.y * c - relative.z * s;
            rotated.z = relative.y * s + relative.z * c;
            break;
        case 'y': // rotate around Y-axis
            rotated.x = relative.x * c + relative.z * s;
            rotated.z = -relative.x * s + relative.z * c;
            break;
        case 'z': // rotate around Z-axis
            rotated.x = relative.x * c - relative.y * s;
            rotated.y = relative.x * s + relative.y * c;
            break;
    }

    // Update position
    position = target + rotated;

    // If keepFacing is true, make camera look at target
    if (keepFacing) {
        Vector3d dir = (target - position).normalized();

        // Extract yaw and pitch from direction
        rotation.y = std::atan2(dir.x, dir.z);       // yaw
        rotation.x = std::asin(-dir.y);              // pitch
        // leave rotation.z (roll) unchanged
    }
}



};



std::vector<Triangle> loadTRIS(const std::string &filename)
{
    std::vector<Triangle> triangles;
    std::ifstream file(filename);

    if (!file)
    {
        std::cerr << "Failed to open file: " << filename << "\n";
        return triangles;
    }

    int count = 0;
    file >> count;

    if (count <= 0)
    {
        std::cerr << "Invalid or missing triangle count.\n";
        return triangles;
    }

    for (int i = 0; i < count; ++i)
    {
        Triangle t;
        if (!(file >> t.v1.x >> t.v1.y >> t.v1.z >> t.v2.x >> t.v2.y >> t.v2.z >> t.v3.x >> t.v3.y >> t.v3.z))
        {
            std::cerr << "Failed to read triangle #" << i << "\n";
            break;
        }
        triangles.push_back(t);
    }

    return triangles;
}

std::vector<Triangle> loadOBJ(const std::string& filename) {
    std::vector<Triangle> triangles;


    //probably objl can be modified to use local types to remove some overhead
    objl::Loader loader;
    std::ifstream file(filename);


    bool isLoaded = loader.LoadFile(filename);
    if(!isLoaded){
        SDL_Log("failed to load OBJ file");
        return triangles;
    }


    for(objl::Mesh mesh: loader.LoadedMeshes){
        for (int j = 0; j < mesh.Indices.size(); j += 3)
			{
				int i1, i2, i3;
                i1 = mesh.Indices[j];
                i2 = mesh.Indices[j + 1];
                i3 = mesh.Indices[j + 2];

                Vector3d v1, v2, v3;

                v1.x = mesh.Vertices[i1].Position.X;
                v1.y = mesh.Vertices[i1].Position.Y;
                v1.z = mesh.Vertices[i1].Position.Z;

                v2.x = mesh.Vertices[i2].Position.X;
                v2.y = mesh.Vertices[i2].Position.Y;
                v2.z = mesh.Vertices[i2].Position.Z;

                v3.x = mesh.Vertices[i3].Position.X;
                v3.y = mesh.Vertices[i3].Position.Y;
                v3.z = mesh.Vertices[i3].Position.Z;

                Triangle tri;
                tri.v1 = v1;
                tri.v2 = v2;
                tri.v3 = v3;

                tri.vn1.x = mesh.Vertices[i1].Normal.X;
                tri.vn1.y = mesh.Vertices[i1].Normal.Y;
                tri.vn1.z = mesh.Vertices[i1].Normal.Z;

                tri.vn2.x = mesh.Vertices[i2].Normal.X;
                tri.vn2.y = mesh.Vertices[i2].Normal.Y;
                tri.vn2.z = mesh.Vertices[i2].Normal.Z;

                tri.vn3.x = mesh.Vertices[i3].Normal.X;
                tri.vn3.y = mesh.Vertices[i3].Normal.Y;
                tri.vn3.z = mesh.Vertices[i3].Normal.Z;

                triangles.push_back(tri);
			}
    }


    return triangles;
}


// Rotate a point around a pivot using Euler angles (rotation.x, rotation.y, rotation.z)
Vector3d rotatePoint(const Vector3d& point, const Vector3d& pivot, const Vector3d& rotation) {
    // Translate point relative to pivot
    float px = point.x - pivot.x;
    float py = point.y - pivot.y;
    float pz = point.z - pivot.z;

    // Rotation angles (radians)
    float cx = cos(rotation.x), sx = sin(rotation.x);
    float cy = cos(rotation.y), sy = sin(rotation.y);
    float cz = cos(rotation.z), sz = sin(rotation.z);

    // Rotation matrix (XYZ order: rotate around X, then Y, then Z)
    float m00 = cy * cz;
    float m01 = -cy * sz;
    float m02 = sy;

    float m10 = sx * sy * cz + cx * sz;
    float m11 = -sx * sy * sz + cx * cz;
    float m12 = -sx * cy;

    float m20 = -cx * sy * cz + sx * sz;
    float m21 = cx * sy * sz + sx * cz;
    float m22 = cx * cy;

    // Apply rotation
    float rx = m00 * px + m01 * py + m02 * pz;
    float ry = m10 * px + m11 * py + m12 * pz;
    float rz = m20 * px + m21 * py + m22 * pz;

    // Translate back to pivot
    return Vector3d(rx + pivot.x, ry + pivot.y, rz + pivot.z);
}

//this bit is ONLY TO BE USED FOR SHOWCASING
//this is not even a believeable light source
struct LambertLight {
    Vector3d direction = Vector3d(-1, 0, 1).normalized(); // should be normalized
    float intensity = 1.0;    // 0..1
};


LambertLight globalLambertLight;


//unholy beast
// Comparator for std::sort
//perhaps calculating all 3 axis for orthographic is causing visual artifacts?
struct TriangleDistanceComparator {
    Vector3d point; // usually the camera position

    TriangleDistanceComparator(const Vector3d& p) : point(p) {}

    bool operator()(Triangle& t1, Triangle& t2) const {
        Vector3d c1 = t1.centroid();
        Vector3d c2 = t2.centroid();

        float d1 = c1.squaredDistance(point);

        float d2 = c2.squaredDistance(point);
        return d1 > d2; 
    }
};


namespace Shader
{
    Triangle* tri;
    uint32_t solidColor, v1color, v2color, v3color;
    float solidintensity, v1intensity, v2intensity, v3intensity;
    int x, y;

    LambertLight l_light;

    uint32_t (*colorFunction)(int, int) = NULL;
    uint32_t (*lightFunction)(int, int) = NULL; 

    void init(){
        tri = NULL;

        solidColor = 0;
        v1color = 0;
        v2color = 0;
        v3color = 0;

        solidintensity = 0;
        v1intensity = 0;
        v2intensity = 0;
        v3intensity = 0;

        colorFunction = NULL;
        lightFunction = NULL;
        
    }

    uint32_t getFlatColor(float alpha, float beta, float gamma){
        return solidColor;
    }

    uint32_t getSmoothColor(float alpha, float beta, float gamma){
        float intensity = alpha * v1intensity + beta * v2intensity + gamma * v3intensity;
        if(intensity > 1.0){ intensity = 1.0; }
        if(intensity < 0.1){ intensity = 0.1; }

        uint8_t r = (solidColor >> 16) & 0xFF;
        uint8_t g = (solidColor >> 8) & 0xFF;
        uint8_t b = solidColor & 0xFF;

        // Apply shading
        r = static_cast<uint8_t>(r * intensity);
        g = static_cast<uint8_t>(g * intensity);
        b = static_cast<uint8_t>(b * intensity);
        return packARGB(255, r, g, b);
    }

    float computeLambertLighting(Vector3d normal) {
        Vector3d n = normal;     // Surface normal
        Vector3d l = l_light.direction;

        float dotNL = n.dot(l);
        float diffuse = std::max(0.0f, dotNL) * l_light.intensity;

        // Add some ambient light so faces in shadow aren't pure black
        float ambient = 0.2f;
        return std::min(1.0f, ambient + diffuse);
    }

    void flatColorRun(){
        int a, r, g, b;
        if(colorMode == ColorMode::RANDOM_COLOR){
            a = tri->a;
            r = tri->r;
            g = tri->g;
            b = tri->b;
        }
        else if(colorMode == ColorMode::SOLID_WHITE){
            a = 255;
            r = 255;
            g = 255;
            b = 255;
        }
        if(lightMode == LightMode::SIMPLE_LAMBERT){
            float intensity = Shader::computeLambertLighting(tri->Normal());
            r *= intensity;
            g *= intensity;
            b *= intensity;
        }
        Shader::solidColor = packARGB(a, r, g, b);
    }

    void smoothColorRun(){
        if(lightMode == LightMode::CONSTANT){
            flatColorRun();
            return;
        }

        int a, r, g, b;
        float i1, i2, i3;
        if(colorMode == ColorMode::RANDOM_COLOR){
            a = tri->a;
            r = tri->r;
            g = tri->g;
            b = tri->b;
        }
        else if(colorMode == ColorMode::SOLID_WHITE){
            a = 255;
            r = 255;
            g = 255;
            b = 255;
        }
        if(lightMode == LightMode::SIMPLE_LAMBERT){
            v1intensity = Shader::computeLambertLighting(tri->vn1);
            v2intensity = Shader::computeLambertLighting(tri->vn2);
            v3intensity = Shader::computeLambertLighting(tri->vn3);
        }
        Shader::solidColor = packARGB(a, r, g, b);

    }

    void colorRun(){
        if(shadeMode == ShadeMode::FLAT) flatColorRun();
        else smoothColorRun();
    }

}

ScreenPoint orthos(const Vector3d &point, const Camera &camera, const ScreenBuffer &screenBuffer) {
    ScreenPoint screenPoint;

    Vector3d relative = point - camera.position;
    Vector3d rotatedPoint = rotatePoint(
        relative, 
        Vector3d(0,0,0), 
        Vector3d(-camera.rotation.x, -camera.rotation.y, -camera.rotation.z)
    );

    screenPoint.x = static_cast<int>(rotatedPoint.x * camera.zoom * camera.aspect + screenBuffer.halfwidth); 
    screenPoint.y = static_cast<int>(-rotatedPoint.y * camera.zoom + screenBuffer.halfheight);
    screenPoint.z = 1/((rotatedPoint.z - camera.near)/(camera.far - camera.near)); 

    return screenPoint;
}

ScreenPoint pers(const Vector3d &point, const Camera &camera, const ScreenBuffer &screenBuffer) {
    ScreenPoint screenPoint;

    Vector3d relative = point - camera.position;
    Vector3d rotatedPoint = rotatePoint(
        relative, 
        Vector3d(0,0,0), 
        Vector3d(-camera.rotation.x, -camera.rotation.y, -camera.rotation.z)
    );

    if (rotatedPoint.z <= 0.001f) 
        rotatedPoint.z = 0.001f;

    float scale = (camera.fov / rotatedPoint.z);

    screenPoint.x = static_cast<int>(rotatedPoint.x * scale * camera.aspect * screenBuffer.halfwidth + screenBuffer.halfwidth);
    screenPoint.y = static_cast<int>(-rotatedPoint.y * scale * screenBuffer.halfheight + screenBuffer.halfheight);
    screenPoint.z = 1/((rotatedPoint.z - camera.near)/(camera.far - camera.near)); 

    return screenPoint;
}


void renderFrame(ScreenBuffer &screenBuffer, Camera &camera, const std::vector<Triangle> &tris){


    //clear the screen
    screenBuffer.fill(0);
    screenBuffer.drainZBuffer();
    

    if(lightMode == LightMode::SIMPLE_LAMBERT) {
        globalLambertLight.direction = camera.getDirection().rotateByAngle(Vector3d(0, 1, 0), 2.0).normalized();
    }


    //projection function
    ScreenPoint (*project)(const Vector3d &point, const Camera &camera, const ScreenBuffer &screenBuffer);
    if(projectionMode == ProjectionMode::PROJ_ORTHOGRAPHIC) project = orthos;
    if(projectionMode == ProjectionMode::PROJ_PERSPECTIVE) project = pers;

    //triangle draw function (no zbuffer)
    void (*drawTriangle)(int ax, int ay, int bx, int by, int cx, int cy, ScreenBuffer &buffer, uint32_t (*color)(float, float, float));
    if(rasterMode == TriangleRasterMethod::SCANLINE) drawTriangle = triangleRasterScanLine; 
    if(rasterMode == TriangleRasterMethod::BOUNDING_BOX) drawTriangle = &triangleRasterBoundingBox;

    //triangle draw function (with zbuffer)
    void (*drawTriangleZBuffer)(int ax, int ay, int bx, int by, int cx, int cy, float az, float bz, float cz, ScreenBuffer &buffer, uint32_t (*color)(float, float, float));
    if(rasterMode == TriangleRasterMethod::SCANLINE) drawTriangleZBuffer = triangleRasterScanline_ZBuffered;
    if(rasterMode == TriangleRasterMethod::BOUNDING_BOX) drawTriangleZBuffer = triangleRasterBoundingBox_ZBuffered;

    //color function
    uint32_t (*colorFunction)(float, float, float);
    if(colorMode != ColorMode::TEXTURE){
        if(shadeMode == ShadeMode::FLAT) colorFunction = Shader::getFlatColor;
        if(shadeMode == ShadeMode::SMOOTH) colorFunction = Shader::getSmoothColor;
    } 


    int draw_count, all_count;
        
    draw_count = 0;
    all_count = 0;

    Vector3d camLoc = camera.position;
    float intensity;


        if( projectionMode == ProjectionMode::PROJ_ORTHOGRAPHIC){

            //not using origin point might cause problem
            camLoc = camLoc + camera.getDirection() * -orthographicHurlUnit;
        }

    std::vector<Triangle> modified_tris;
    if(renderMode == RenderMode::RENDER_PAINTERS){
        modified_tris = tris;

        std::sort(modified_tris.begin(), modified_tris.end(), TriangleDistanceComparator(camLoc));
    }
    else{
        modified_tris = std::move(tris); 
    }

    for(Triangle t: modified_tris){
        all_count++;

        if(backfaceCulling){
            if(camera.getDirection().dot(t.Normal()) >= 0.1){ continue; }
        }
        
        draw_count++;

        //
        Shader::init();
        Shader::tri = &t;
        Shader::colorRun();
        //


        ScreenPoint p1, p2, p3;
        p1 = project(t.v1, camera, screenBuffer);
        p2 = project(t.v2, camera, screenBuffer);
        p3 = project(t.v3, camera, screenBuffer);
        

        if(renderMode == RenderMode::RENDER_WIREFRAME){
            drawLine(p1.x, p1.y, p2.x, p2.y, screenBuffer, 0xFFFFFFFF);
            drawLine(p2.x, p2.y, p3.x, p3.y, screenBuffer, 0xFFFFFFFF);
            drawLine(p3.x, p3.y, p1.x, p1.y, screenBuffer, 0xFFFFFFFF);
            continue;
        }

        //solid lights and color


        if(renderMode == RenderMode::RENDER_ZBUFFER_MODEL || renderMode == RenderMode::RENDER_ZBUFFER_BUFFER){
            drawTriangleZBuffer(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p1.z, p2.z, p3.z, screenBuffer, colorFunction);
            continue;            
        }        
        if(renderMode == RenderMode::RENDER_PAINTERS || renderMode == RenderMode::RENDER_NO_DEPTH){
            drawTriangle(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, screenBuffer, colorFunction);
        }

    }


    if(renderMode == RenderMode::RENDER_ZBUFFER_BUFFER){
        float minz = *std::min_element(screenBuffer.zBuffer.begin(), screenBuffer.zBuffer.end());
        for(int i = 0; i < screenBuffer.buffer.size(); i++){
            float z = screenBuffer.zBuffer[i];
            if(z > minz + ZDEPTH_VISIBLE){ screenBuffer.buffer[i] = 0; continue; }
            float zNorm = 1 - (z - minz)/ (ZDEPTH_VISIBLE - minz); // I WILL put MAGIC NUMBERS in code
            
            screenBuffer.buffer[i] = packARGB(255, 255 * zNorm, 255 * zNorm, 255 * zNorm);
        }
    }

}

int main() {

    ScreenBuffer screenBuffer(1000, 1000);
    
    std::vector<Triangle> tris;
    std::vector<std::vector<Triangle>> models;


    tris = loadOBJ("african_head.obj");


    models.push_back(tris);

    Camera mainCamera;
    mainCamera.aspect = screenBuffer.width / screenBuffer.height;
    mainCamera.zoom = 250.0f;
    mainCamera.position.y = 0;
    mainCamera.position.z = -2;
    

    SDL_Window* window = SDL_CreateWindow("Hello SDL3",
                                          screenBuffer.width, screenBuffer.height, 0);
    if (!window) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << "\n";
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
    if (!renderer) {
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << "\n";
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    SDL_Texture* framebuffer = SDL_CreateTexture(
        renderer,
        SDL_PIXELFORMAT_ARGB8888,
        SDL_TEXTUREACCESS_STREAMING,
        screenBuffer.width, screenBuffer.height);


    

    bool running = true;
    SDL_Event event;

    while (running) {

        Uint64 start = SDL_GetPerformanceCounter();
        // Handle events
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                running = false;
            }
        }




        
        // Clear screen with a color (red, green, blue, alpha)
        SDL_RenderClear(renderer);

        renderFrame(screenBuffer, mainCamera, tris);


        
        // Present the frame
        SDL_UpdateTexture(framebuffer, nullptr, screenBuffer.buffer.data(), screenBuffer.width * sizeof(uint32_t));

        // render
        SDL_RenderClear(renderer);
        SDL_RenderTexture(renderer, framebuffer, nullptr, nullptr);
        SDL_RenderPresent(renderer);

	    Uint64 end = SDL_GetPerformanceCounter();

	    float delta = (end - start) / (float)SDL_GetPerformanceFrequency();

        mainCamera.orbitAroundPoint(Vector3d(), 0.4f * delta, 'y', true);
        
        //lambert rotation
        Shader::l_light.direction = mainCamera.getDirection().cross(Vector3d(0, -1, 0));
        //Shader::l_light.direction.z -= 0.4;


        std::string title = "fps: " + std::to_string(1.0f / delta);
	    SDL_SetWindowTitle(window , title.c_str());
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
