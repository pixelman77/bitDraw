#include <SDL3/SDL.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

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
    RENDER_ZBUFFER_BUFFER
};


//buffer input gets a seizure with these and I have no idea why
uint32_t SCREEN_WIDTH = 800;
uint32_t SCREEN_HEIGHT = 600;
float HALFWIDTH = SCREEN_WIDTH / 2;
float HALFHEIGHT = SCREEN_HEIGHT / 2;


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

Vector3d cross(const Vector3d &a, const Vector3d &b)
{
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}



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


Vector2di orthos(Vector3d point, Camera &camera, ScreenBuffer &screenBuffer){
    
    Vector2di screenPoint;
    Vector3d relative = point - camera.position;
    Vector3d rotatedPoint = rotatePoint(relative, Vector3d(0,0,0), Vector3d(-camera.rotation.x, -camera.rotation.y, -camera.rotation.z));

    screenPoint.x = rotatedPoint.x * camera.zoom * camera.aspect + screenBuffer.halfwidth; 
    screenPoint.y = -rotatedPoint.y * camera.zoom + screenBuffer.halfheight;

    return screenPoint;
}

int main() {


    ScreenBuffer screenBuffer(800, 600);
    std::vector<Triangle> tris = loadTRIS("utah_teapot.tris");
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
        screenBuffer.fill(0);

        int c = 0;

        for(Triangle t: tris){
            Vector2di p1, p2, p3;
            p1 = orthos(t.v1, mainCamera, screenBuffer);
            p2 = orthos(t.v2, mainCamera, screenBuffer);
            p3 = orthos(t.v3, mainCamera, screenBuffer);

            screenBuffer.line(p1.x, p1.y, p2.x, p2.y, 0xFFFFFFFF);
            screenBuffer.line(p2.x, p2.y, p3.x, p3.y, 0xFFFFFFFF);
            screenBuffer.line(p3.x, p3.y, p1.x, p1.y, 0xFFFFFFFF);

            SDL_Log(std::to_string(c).c_str());
            c++;
        }
        


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
