#pragma once
#include <cstdint>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstring>

#pragma pack(push,1)
struct TGAHeader {
    std::uint8_t  idlength = 0;
    std::uint8_t  colormaptype = 0;
    std::uint8_t  datatypecode = 0;
    std::uint16_t colormaporigin = 0;
    std::uint16_t colormaplength = 0;
    std::uint8_t  colormapdepth = 0;
    std::uint16_t x_origin = 0;
    std::uint16_t y_origin = 0;
    std::uint16_t width = 0;
    std::uint16_t height = 0;
    std::uint8_t  bitsperpixel = 0;
    std::uint8_t  imagedescriptor = 0;
};
#pragma pack(pop)


struct TGAColor {
    std::uint8_t bgra[4] = {0,0,0,0};
    std::uint8_t bytespp = 4;

    TGAColor() = default;

    TGAColor(std::uint8_t R, std::uint8_t G, std::uint8_t B, std::uint8_t A = 255) {
        bgra[0] = B; bgra[1] = G; bgra[2] = R; bgra[3] = A;
    }

    std::uint8_t& operator[](const int i) { return bgra[i]; }
};


struct TGAImage {
    enum Format { GRAYSCALE = 1, RGB = 3, RGBA = 4 };

    TGAImage() = default;

    TGAImage(int w, int h, int bpp)
        : w(w), h(h), bpp(bpp), data(w*h*bpp, 0) {}

    TGAImage(int w, int h, int bpp, TGAColor c)
        : w(w), h(h), bpp(bpp), data(w*h*bpp, 0)
    {
        for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
                set(x, y, c);
    }

    // ---------------- READ TGA ----------------
    bool read_tga_file(const std::string filename) {
        std::ifstream in(filename, std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "can't open file " << filename << "\n";
            return false;
        }

        TGAHeader header {};
        in.read(reinterpret_cast<char*>(&header), sizeof(header));
        if (!in.good()) {
            std::cerr << "error reading header\n";
            return false;
        }

        w = header.width;
        h = header.height;
        bpp = header.bitsperpixel >> 3;

        if (w <= 0 || h <= 0 || (bpp != GRAYSCALE && bpp != RGB && bpp != RGBA)) {
            std::cerr << "bad bpp or dimension\n";
            return false;
        }

        size_t nbytes = size_t(w) * h * bpp;
        data.assign(nbytes, 0);

        // Uncompressed
        if (header.datatypecode == 2 || header.datatypecode == 3) {
            in.read(reinterpret_cast<char*>(data.data()), nbytes);
            if (!in.good()) return false;
        }
        // Compressed RLE
        else if (header.datatypecode == 10 || header.datatypecode == 11) {
            if (!load_rle_data(in)) return false;
        }
        else {
            std::cerr << "unknown TGA type: " << int(header.datatypecode) << "\n";
            return false;
        }

        // Orientation
        if (!(header.imagedescriptor & 0x20))
            flip_vertically();
        if (header.imagedescriptor & 0x10)
            flip_horizontally();

        return true;
    }


    // --------------- WRITE TGA ----------------
    bool write_tga_file(const std::string filename, bool vflip = true, bool rle = true) const {
        constexpr std::uint8_t developer_area_ref[4] = {0,0,0,0};
        constexpr std::uint8_t extension_area_ref[4] = {0,0,0,0};
        constexpr std::uint8_t footer[18] = {
            'T','R','U','E','V','I','S','I','O','N','-','X','F','I','L','E','.','\0'
        };

        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) {
            std::cerr << "can't open file " << filename << "\n";
            return false;
        }

        TGAHeader header {};
        header.width = w;
        header.height = h;
        header.bitsperpixel = bpp * 8;
        header.datatypecode = (bpp == 1 ? (rle ? 11 : 3)
                                        : (rle ? 10 : 2));
        header.imagedescriptor = vflip ? 0x00 : 0x20;

        if (!out.write(reinterpret_cast<const char*>(&header), sizeof(header)))
            return false;

        if (!rle) {
            if (!out.write(reinterpret_cast<const char*>(data.data()), size_t(w)*h*bpp))
                return false;
        }
        else {
            if (!unload_rle_data(out))
                return false;
        }

        if (!out.write(reinterpret_cast<const char*>(developer_area_ref), sizeof(developer_area_ref)))
            return false;
        if (!out.write(reinterpret_cast<const char*>(extension_area_ref), sizeof(extension_area_ref)))
            return false;
        if (!out.write(reinterpret_cast<const char*>(footer), sizeof(footer)))
            return false;

        return true;
    }


    // --------------- PIXEL ACCESS ----------------
    TGAColor get(int x, int y) const {
        if (!data.size() || x < 0 || y < 0 || x >= w || y >= h)
            return {};

        TGAColor c;
        c.bytespp = bpp;

        const uint8_t* p = data.data() + (x + y * w) * bpp;
        for (int i = 0; i < bpp; i++)
            c.bgra[i] = p[i];

        return c;
    }

    void set(int x, int y, const TGAColor& c) {
        if (!data.size() || x < 0 || y < 0 || x >= w || y >= h)
            return;

        memcpy(data.data() + (x + y * w) * bpp, c.bgra, bpp);
    }


    // --------------- FLIP IMAGE ----------------
    void flip_horizontally() {
        for (int x = 0; x < w / 2; x++)
            for (int y = 0; y < h; y++)
                for (int b = 0; b < bpp; b++)
                    std::swap(
                        data[(x + y*w)*bpp + b],
                        data[((w - 1 - x) + y*w)*bpp + b]
                    );
    }

    void flip_vertically() {
        for (int y = 0; y < h / 2; y++)
            for (int x = 0; x < w; x++)
                for (int b = 0; b < bpp; b++)
                    std::swap(
                        data[(x + y*w)*bpp + b],
                        data[(x + (h - 1 - y)*w)*bpp + b]
                    );
    }


    // --------------- INFO ----------------
    int width() const { return w; }
    int height() const { return h; }


private:

    // --------- LOAD RLE ----------
    bool load_rle_data(std::ifstream& in) {
        size_t pixelcount = size_t(w) * h;
        size_t currentpixel = 0;
        size_t currentbyte = 0;

        TGAColor colorbuffer;

        while (currentpixel < pixelcount) {
            uint8_t chunk = in.get();
            if (!in.good()) return false;

            if (chunk < 128) {
                chunk++;
                for (int i = 0; i < chunk; i++) {
                    in.read(reinterpret_cast<char*>(colorbuffer.bgra), bpp);
                    if (!in.good()) return false;

                    for (int t = 0; t < bpp; t++)
                        data[currentbyte++] = colorbuffer.bgra[t];

                    if (++currentpixel > pixelcount) return false;
                }
            }
            else {
                chunk -= 127;
                in.read(reinterpret_cast<char*>(colorbuffer.bgra), bpp);
                if (!in.good()) return false;

                for (int i = 0; i < chunk; i++) {
                    for (int t = 0; t < bpp; t++)
                        data[currentbyte++] = colorbuffer.bgra[t];

                    if (++currentpixel > pixelcount) return false;
                }
            }
        }
        return true;
    }


    // --------- WRITE RLE ----------
    bool unload_rle_data(std::ofstream& out) const {
        const uint8_t max_chunk = 128;
        size_t npixels = size_t(w) * h;
        size_t curpix = 0;

        while (curpix < npixels) {
            size_t chunkstart = curpix * bpp;
            size_t curbyte = curpix * bpp;

            uint8_t run_length = 1;
            bool raw = true;

            while (curpix + run_length < npixels && run_length < max_chunk) {
                bool equal = true;
                for (int t = 0; t < bpp; t++)
                    equal = equal && (data[curbyte + t] == data[curbyte + t + bpp]);

                curbyte += bpp;

                if (run_length == 1)
                    raw = !equal;

                if (raw && equal) {
                    run_length--;
                    break;
                }
                if (!raw && !equal)
                    break;

                run_length++;
            }

            curpix += run_length;

            out.put(raw ? run_length - 1 : run_length + 127);
            if (!out.good()) return false;

            out.write(reinterpret_cast<const char*>(data.data() + chunkstart),
                      raw ? run_length * bpp : bpp);

            if (!out.good()) return false;
        }

        return true;
    }


    int w = 0, h = 0;
    std::uint8_t bpp = 0;
    std::vector<std::uint8_t> data;
};

