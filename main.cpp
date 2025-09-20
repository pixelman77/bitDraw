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

int SCREEN_WIDTH = 800;
int SCREEN_HEIGHT = 600;
float HALFWIDTH = SCREEN_WIDTH / 2;
float HALFHEIGHT = SCREEN_HEIGHT / 2;

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


Vector2di orthos(Vector3d point, Camera &camera){
    
    Vector3d rotatedPoint = rotatePoint(point, camera.position, camera.rotation);
    Vector2di screenPoint;

    screenPoint.x = rotatedPoint.x * camera.zoom * camera.aspect + HALFWIDTH; 
    screenPoint.y = -rotatedPoint.y * camera.zoom + HALFHEIGHT;

    return screenPoint;
}

int main() {

    std::vector<Triangle> tris = loadTRIS("utah_teapot.tris");
    Camera mainCamera;

    SDL_Window* window = SDL_CreateWindow("Hello SDL3",
                                          SCREEN_WIDTH, SCREEN_HEIGHT,
                                          SDL_WINDOW_VULKAN);
    if (!window) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << "\n";
        SDL_Quit();
        return 1;
    }

    // Create a renderer so we can draw things
    SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
    if (!renderer) {
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << "\n";
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    //little testing bit 
    for(Triangle t : tris){
        t.v1.z += 5.0f;
        t.v2.z += 5.0f;
        t.v3.z += 5.0f;
    }
    Vector3d pivot = tris[10].v2;

    //end of little testing bit


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
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); 
        SDL_RenderClear(renderer);

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        for(Triangle t: tris){
            Vector2di p1, p2, p3;
            p1 = orthos(t.v1, mainCamera);
            p2 = orthos(t.v2, mainCamera);
            p3 = orthos(t.v3, mainCamera);

            SDL_RenderLine(renderer, p1.x, p1.y, p2.x, p2.y);
            SDL_RenderLine(renderer, p2.x, p2.y, p3.x, p3.y);
            SDL_RenderLine(renderer, p3.x, p3.y, p1.x, p1.y);

        }
        

        mainCamera.rotation.y += 0.005f;
        // Present the frame
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
