//
// Created by minaj on 10/4/2025.
//

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "structs/vertex.h"
#include "geometries/geometry.h"

class Triangle : public Geometry {
public:
    Triangle(const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr);
private:
    void FillVertexBuffer(const std::vector<Vertex>& vertices);
    void FillIndexBuffer(const std::vector<Vertex>& vertices);

    float* m_vertexBuffer = nullptr;
    unsigned* m_indexBuffer = nullptr;
};



#endif //TRIANGLE_H
