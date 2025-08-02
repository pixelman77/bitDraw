#include <SFML/Graphics.hpp>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>


float moveSpeed = 0.1f;
float rotSpeed = 0.08f;


struct Vector2d {
    float x, y;

    // arithmetc logic?
};

struct IVector2d{
    int x, y;
};

struct Vector3d {
    float x, y, z;

    Vector3d(float x = 0, float y = 0, float z = 0){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector3d operator+(const Vector3d& rhs) const {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    Vector3d operator-(const Vector3d& rhs) const {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }
};

struct Triangle {
    Vector3d v1, v2, v3;

    Triangle(){
        v1 = Vector3d();
        v2 = Vector3d();
        v3 = Vector3d();
    }

};

struct Camera {
    Vector3d position, rotation;
    float fov;    // Field of view
    float aspect; // Aspect ratio (width / height)

    Camera(){
        position = Vector3d();
        rotation = Vector3d();
        fov = 1.0f;
        aspect = 800.0f / 600.0f;
    }
};

void moveCamera(Camera& camera){

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::A)) {
            camera.position.x -= cos(camera.rotation.y) * moveSpeed;
            camera.position.z -= sin(camera.rotation.y) * moveSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::D)) {
            camera.position.x += cos(camera.rotation.y) * moveSpeed;
            camera.position.z += sin(camera.rotation.y) * moveSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::W)) {
            camera.position.x -= sin(camera.rotation.y) * moveSpeed;
            camera.position.z += cos(camera.rotation.y) * moveSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::S)) {
            camera.position.x += sin(camera.rotation.y) * moveSpeed;
            camera.position.z -= cos(camera.rotation.y) * moveSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q)) {
            camera.position.y += moveSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::E)) {
            camera.position.y -= moveSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) {
            camera.rotation.y += rotSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) {
            camera.rotation.y -= rotSpeed;
        }

        //rotation is not local
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) {
            camera.rotation.x += rotSpeed;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) {
            camera.rotation.x -= rotSpeed;
        }
}

sf::Vector2f perspectiveProject(const Vector3d& point, const Camera& cam) {
    // Translate to camera space
    Vector3d relative = point - cam.position;

    // Apply pitch (X-axis)
    float cosPitch = std::cos(cam.rotation.x);
    float sinPitch = std::sin(cam.rotation.x);
    float y1 = relative.y * cosPitch - relative.z * sinPitch;
    float z1 = relative.y * sinPitch + relative.z * cosPitch;
    float x1 = relative.x;

    // Apply yaw (Y-axis)
    float cosYaw = std::cos(cam.rotation.y);
    float sinYaw = std::sin(cam.rotation.y);
    float x2 = x1 * cosYaw + z1 * sinYaw;
    float z2 = -x1 * sinYaw + z1 * cosYaw;
    float y2 = y1;

    // Apply roll (Z-axis)
    float cosRoll = std::cos(cam.rotation.z);
    float sinRoll = std::sin(cam.rotation.z);
    float x3 = x2 * cosRoll - y2 * sinRoll;
    float y3 = x2 * sinRoll + y2 * cosRoll;

    // Perspective projection
    if (z2 <= 0.001f) z2 = 0.001f; // avoid division by zero

    float scale = cam.fov / z2;
    float screenX = x3 * scale * cam.aspect * 300.0f + 400.0f; // assuming 800x600 screen
    float screenY = -y3 * scale * 300.0f + 300.0f;

    return {screenX, screenY};
}

std::vector<Triangle> loadTRIS(const std::string& filename) {
    std::vector<Triangle> triangles;
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return triangles;
    }

    int count = 0;
    file >> count;

    if (count <= 0) {
        std::cerr << "Invalid or missing triangle count.\n";
        return triangles;
    }

    for (int i = 0; i < count; ++i) {
        Triangle t;
        if (!(file >> t.v1.x >> t.v1.y >> t.v1.z
                  >> t.v2.x >> t.v2.y >> t.v2.z
                  >> t.v3.x >> t.v3.y >> t.v3.z)) {
            std::cerr << "Failed to read triangle #" << i << "\n";
            break;
        }
        triangles.push_back(t);
    }

    return triangles;
}


int main() {    
    float speed = 0.5f;
    Camera stCamera;

    stCamera.position = Vector3d(0, 0, 0);
    //stCamera.rotation = Vector3d(-0.6, -2.3, 0);

    std::vector<Triangle> triangles = loadTRIS("model.tris");

    sf::RenderWindow window(sf::VideoMode(800, 600), "TRIS Viewer (SFML Wireframe)");
    window.setFramerateLimit(60);

    while (window.isOpen()) {
        sf::Event e;
        while (window.pollEvent(e)) {
            if (e.type == sf::Event::Closed)
                window.close();
        }

        moveCamera(stCamera);
        window.clear(sf::Color::Black);

        for (const Triangle& tri : triangles) {

            sf::Vector2f p1 = perspectiveProject(tri.v1, stCamera);
            sf::Vector2f p2 = perspectiveProject(tri.v2, stCamera);
            sf::Vector2f p3 = perspectiveProject(tri.v3, stCamera);

            sf::Vertex lines[] = {
                sf::Vertex(p1, sf::Color::White),
                sf::Vertex(p2, sf::Color::White),
                sf::Vertex(p3, sf::Color::White),
                sf::Vertex(p1, sf::Color::White)
            };

            window.draw(lines, 4, sf::LineStrip);
        }

        window.display();

    }

}
