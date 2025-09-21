#include <SDL3/SDL.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <bits/stdc++.h>

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

//default options
ProjectionMode projectionMode = ProjectionMode::PROJ_PERSPECTIVE;
RenderMode renderMode = RenderMode::RENDER_PAINTERS;
TriangleRasterMethod rasterMode = TriangleRasterMethod::SCANLINE;
bool backfaceCulling = true;


//buffer input gets a seizure with these and I have no idea why
uint32_t SCREEN_WIDTH = 800;
uint32_t SCREEN_HEIGHT = 600;
float HALFWIDTH = SCREEN_WIDTH / 2;
float HALFHEIGHT = SCREEN_HEIGHT / 2;

uint32_t packARGB(uint8_t a, uint8_t r, uint8_t g, uint8_t b) {
    return (static_cast<uint32_t>(a) << 24) |
           (static_cast<uint32_t>(r) << 16) |
           (static_cast<uint32_t>(g) << 8)  |
           (static_cast<uint32_t>(b));
}

double signed_triangle_area(int ax, int ay, int bx, int by, int cx, int cy) {
    return .5*((by-ay)*(bx+ax) + (cy-by)*(cx+bx) + (ay-cy)*(ax+cx));
}

struct ScreenBuffer
{
    // note to self: w of h/x of y
    // down growing
    std::vector<uint32_t> buffer;
    int width, height;
    float halfwidth, halfheight;

    ScreenBuffer(uint32_t iwidth, uint32_t iheigth){

        this->width = iwidth;
        this->height = iheigth;
        halfheight = height / 2;
        halfwidth = width / 2;
        buffer = std::vector<uint32_t>(width * height, 0);
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

    void line( int x0, int y0, int x1, int y1, uint32_t color) {
        int dx = std::abs(x1 - x0);
        int dy = std::abs(y1 - y0);
        int sx = (x0 < x1) ? 1 : -1;
        int sy = (y0 < y1) ? 1 : -1;
        int err = dx - dy;

        while (true) {
            set(x0, y0, color); // draw the pixel

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


    void triangleRasterScanLine(int ax, int ay, int bx, int by, int cx, int cy, u_int32_t color) {
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
                    set(x, y, color);
            }
        }
        if (by != cy) { // if the upper half is not degenerate
            int segment_height = cy - by;
            for (int y=by; y<=cy; y++) { // sweep the horizontal line from by to cy
                int x1 = ax + ((cx - ax)*(y - ay)) / total_height;
                int x2 = bx + ((cx - bx)*(y - by)) / segment_height;
                for (int x=std::min(x1,x2); x<std::max(x1,x2); x++)  // draw a horizontal line
                    set(x, y, color);
            }
        }
    }

    void triangleRasterBoundingBox(int ax, int ay, int bx, int by, int cx, int cy, uint32_t color) {
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
                set(x, y, color);
            }
        }
    }

};

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

struct Vector2di
{
    int x, y;
    Vector2di(int x = 0, int y = 0)
    {
        this->x = x;
        this->y = y;
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

    Vector3d operator+(const Vector3d &rhs) const
    {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    Vector3d operator-(const Vector3d &rhs) const
    {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }

    Vector3d normalized()
    {
        float len = std::sqrt(x * x + y * y + z * z);
        if (len == 0)
            return {0, 0, 0};
        return {x / len, y / len, z / len};
    }
};

struct Triangle
{
    Vector3d v1, v2, v3;

    Triangle()
    {
        v1 = Vector3d();
        v2 = Vector3d();
        v3 = Vector3d();
    }
};

struct Camera
{
    Vector3d position, rotation;
    float fov;    // for perspective
    float zoom;   // for zoom
    float aspect; // Aspect ratio (width / height)

    Camera()
    {
        position = Vector3d();
        rotation = Vector3d();
        fov = 1.0f;
        zoom = 80.0f;
        aspect = SCREEN_WIDTH / SCREEN_HEIGHT;
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

float squaredDistance(const Vector3d& a, const Vector3d& b) {
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

Vector3d centroid(const Triangle& t) {
    return {
        (t.v1.x + t.v2.x + t.v3.x) / 3.0f,
        (t.v1.y + t.v2.y + t.v3.y) / 3.0f,
        (t.v1.z + t.v2.z + t.v3.z) / 3.0f
    };
}

//unholy beast
// Comparator for std::sort
struct TriangleComparator {
    Vector3d point; // usually the camera position

    TriangleComparator(const Vector3d& p) : point(p) {}

    bool operator()(const Triangle& t1, const Triangle& t2) const {
        Vector3d c1 = {
            (t1.v1.x + t1.v2.x + t1.v3.x) / 3.0f,
            (t1.v1.y + t1.v2.y + t1.v3.y) / 3.0f,
            (t1.v1.z + t1.v2.z + t1.v3.z) / 3.0f
        };
        Vector3d c2 = {
            (t2.v1.x + t2.v2.x + t2.v3.x) / 3.0f,
            (t2.v1.y + t2.v2.y + t2.v3.y) / 3.0f,
            (t2.v1.z + t2.v2.z + t2.v3.z) / 3.0f
        };

        float d1 = (c1.x - point.x)*(c1.x - point.x) +
                   (c1.y - point.y)*(c1.y - point.y) +
                   (c1.z - point.z)*(c1.z - point.z);

        float d2 = (c2.x - point.x)*(c2.x - point.x) +
                   (c2.y - point.y)*(c2.y - point.y) +
                   (c2.z - point.z)*(c2.z - point.z);

        return d1 > d2; 
    }
};


Vector2di orthos(const Vector3d &point, const Camera &camera, const ScreenBuffer &screenBuffer){
    
    Vector2di screenPoint;
    Vector3d relative = point - camera.position;
    Vector3d rotatedPoint = rotatePoint(relative, Vector3d(0,0,0), Vector3d(-camera.rotation.x, -camera.rotation.y, -camera.rotation.z));

    screenPoint.x = rotatedPoint.x * camera.zoom * camera.aspect + screenBuffer.halfwidth; 
    screenPoint.y = -rotatedPoint.y * camera.zoom + screenBuffer.halfheight;

    return screenPoint;
}

void renderFrame(ScreenBuffer &screenBuffer, const Camera &camera, const std::vector<std::vector<Triangle>> &models){


    //clear the screen
    screenBuffer.fill(0);
    
    //projection function
    Vector2di (*projFun)(const Vector3d &point, const Camera &camera, const ScreenBuffer &screenBuffer);
    projFun = orthos;

    //triangle draw function
    void (ScreenBuffer::*triDrFunc)(int ax, int ay, int bx, int by, int cx, int cy, uint32_t color);
    if(rasterMode == TriangleRasterMethod::SCANLINE){ triDrFunc = &ScreenBuffer::triangleRasterScanLine; }
    if(rasterMode == TriangleRasterMethod::BOUNDING_BOX){ triDrFunc = &ScreenBuffer::triangleRasterBoundingBox; }



    int tr_count;
    for(std::vector<Triangle> tris : models){
        
        tr_count = 0;

        std::vector<Triangle> modified_tris;
        if(renderMode == RenderMode::RENDER_PAINTERS){
            modified_tris = tris;

            std::sort(modified_tris.begin(), modified_tris.end(), TriangleComparator(camera.position));
        }
        else{
            modified_tris = std::move(tris);
        }

        for(Triangle t: modified_tris){

            if(backfaceCulling){
                //backface culling logic goes here (no shit)
            }

            Vector2di p1, p2, p3;
            p1 = projFun(t.v1, camera, screenBuffer);
            p2 = projFun(t.v2, camera, screenBuffer);
            p3 = projFun(t.v3, camera, screenBuffer);

            //placeholder code until I decide how to handle loose vertexes
            if(renderMode == RenderMode::RENDER_WIREFRAME){

                screenBuffer.line(p1.x, p1.y, p2.x, p2.y, 0xFFFFFFFF);
                screenBuffer.line(p2.x, p2.y, p3.x, p3.y, 0xFFFFFFFF);
                screenBuffer.line(p3.x, p3.y, p1.x, p1.y, 0xFFFFFFFF);
            }
            if(renderMode == RenderMode::RENDER_NO_DEPTH || renderMode == RenderMode::RENDER_PAINTERS){
                uint8_t r, g, b;
                srand(tr_count);
                r = rand() % 256;
                srand(tr_count + 10);
                g = rand() % 256;
                srand(tr_count + 20);
                b = rand() % 256;
                uint32_t color = packARGB(255, r, g, b);
                (screenBuffer.*triDrFunc)(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, color);
            }


            tr_count++;
        }
    }
        

}

int main() {


    ScreenBuffer screenBuffer(800, 600);
    
    std::vector<Triangle> tris = loadTRIS("utah_teapot.tris");
    std::vector<std::vector<Triangle>> models;

    models.push_back(tris);

    Camera mainCamera;
    mainCamera.aspect = screenBuffer.width / screenBuffer.height;
    mainCamera.zoom = 90.0f;
    mainCamera.position.y = 2;
    

    SDL_Window* window = SDL_CreateWindow("Hello SDL3",
                                          SCREEN_WIDTH, SCREEN_HEIGHT, 0);
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

        renderFrame(screenBuffer, mainCamera, models);


        mainCamera.rotation.y += 0.005f;
        // Present the frame
        SDL_UpdateTexture(framebuffer, nullptr, screenBuffer.buffer.data(), screenBuffer.width * sizeof(uint32_t));

        // render
        SDL_RenderClear(renderer);
        SDL_RenderTexture(renderer, framebuffer, nullptr, nullptr);
        SDL_RenderPresent(renderer);

	    Uint64 end = SDL_GetPerformanceCounter();

	    float elapsed = (end - start) / (float)SDL_GetPerformanceFrequency();
        std::string title = "fps: " + std::to_string(1.0f / elapsed);
	    SDL_SetWindowTitle(window , title.c_str());
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
