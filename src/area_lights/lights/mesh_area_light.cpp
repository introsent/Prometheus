#include "mesh_area_light.h"
#include "mesh.h"
#include "math_helpers.h"
#include <algorithm>
#include <iostream>
#include "acceleration/bvh_builder.h"
#include "render/scene_manager.h"
#include "samplers/visibility_aware_sampler.h"

// static member initialization
std::unique_ptr<SpatialVisibilityCache> MeshAreaLight::s_visibilityCache = nullptr;

MeshAreaLight::MeshAreaLight(
    unsigned int meshIndex,
    const glm::vec3& emission,
    float intensity,
    SceneManager* scene)
    : m_meshIndex(meshIndex)
    , m_scene(scene)
    , m_emission(emission)
    , m_intensity(intensity)
    , m_strategy(SamplingStrategy::Uniform)
    , m_totalArea(0.0f)
    , m_rootNodeIndex(-1)
    , m_lastSampledNode(-1)
{
    extractTrianglesFromMesh();
    buildAreaCDF();

    std::cout << "MeshAreaLight created with " << m_triangles.size()
              << " triangles, total area: " << m_totalArea << "\n";
}

MeshAreaLight::~MeshAreaLight() = default;


/// Initialization
void MeshAreaLight::extractTrianglesFromMesh() {
    Mesh* mesh = m_scene->getMesh(m_meshIndex);
    if (!mesh) {
        std::cerr << "ERROR: Mesh not found: " << m_meshIndex << "\n";
        return;
    }

    const auto& vertices = mesh->getOriginalVertices();
    const auto& indices = mesh->getIndices();

    const size_t numTriangles = indices.size() / 3;
    m_triangles.reserve(numTriangles);

    for (size_t i = 0; i < indices.size(); i += 3) {
        // get vertex indices
        const uint32_t idx0 = indices[i];
        const uint32_t idx1 = indices[i + 1];
        const uint32_t idx2 = indices[i + 2];

        // extract positions
        const glm::vec3 v0(vertices[idx0].position.x,
                          vertices[idx0].position.y,
                          vertices[idx0].position.z);
        const glm::vec3 v1(vertices[idx1].position.x,
                          vertices[idx1].position.y,
                          vertices[idx1].position.z);
        const glm::vec3 v2(vertices[idx2].position.x,
                          vertices[idx2].position.y,
                          vertices[idx2].position.z);

        // compute triangle properties
        const glm::vec3 edge1 = v1 - v0;
        const glm::vec3 edge2 = v2 - v0;
        const glm::vec3 crossProd = glm::cross(edge1, edge2);
        const float area = 0.5f * glm::length(crossProd);

        // skip degenerate triangles
        if (area <= 1e-8f) {
            continue;
        }

        const glm::vec3 normal = glm::normalize(crossProd);

        // create triangle data
        TriangleData tri;
        tri.v0 = v0;
        tri.v1 = v1;
        tri.v2 = v2;
        tri.normal = normal;
        tri.area = area;
        tri.intensity = m_intensity;

        // create samplers
        tri.uniformSampler = std::make_unique<UniformTriangleSampler>(
            v0, v1, v2, normal, area, m_emission, m_intensity);

        tri.areaImportanceSampler = std::make_unique<AreaImportanceTriangleSampler>(
            v0, v1, v2, normal, area, m_emission, m_intensity);

        m_triangles.push_back(std::move(tri));
        m_totalArea += area;
    }
}

void MeshAreaLight::buildAreaCDF() {
    if (m_triangles.empty() || m_totalArea <= 0.0f) {
        return;
    }

    m_areaCDF.reserve(m_triangles.size());

    float cumulativeArea = 0.0f;
    for (const auto& tri : m_triangles) {
        cumulativeArea += tri.area;
        m_areaCDF.push_back(cumulativeArea / m_totalArea);
    }

    // ensure exact 1.0 at end (avoid floating point errors)
    m_areaCDF.back() = 1.0f;
}


/// Main sampling interface
AreaLightSample MeshAreaLight::sample(
    const glm::vec3& shadingPoint,
    float u1, float u2, float u3) const {

    switch (m_strategy) {
        case SamplingStrategy::Uniform:
            return sampleUniform(shadingPoint, u1, u2, u3);

        case SamplingStrategy::AreaImportance:
            return sampleAreaImportance(shadingPoint, u1, u2, u3);

        case SamplingStrategy::HierarchicalFlux:
            return sampleHierarchicalFlux(shadingPoint, u1, u2, u3);

        case SamplingStrategy::VisibilityAwareHierarchical:
            return sampleVisibilityAware(shadingPoint, u1, u2, u3);

        default:
            return sampleUniform(shadingPoint, u1, u2, u3);
    }
}

float MeshAreaLight::pdf(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    switch (m_strategy) {
        case SamplingStrategy::Uniform:
            return pdfUniform(shadingPoint, lightPoint);

        case SamplingStrategy::AreaImportance:
            return pdfAreaImportance(shadingPoint, lightPoint);

        case SamplingStrategy::HierarchicalFlux:
            return pdfHierarchicalFlux(shadingPoint, lightPoint);

        case SamplingStrategy::VisibilityAwareHierarchical:
            // visibility-aware pdf is computed by the sampler
            return pdfHierarchicalFlux(shadingPoint, lightPoint);

        default:
            return pdfUniform(shadingPoint, lightPoint);
    }
}


/// Uniform sampling
// select triangle proportional to area, then sample uniformly on triangle

AreaLightSample MeshAreaLight::sampleUniform(
    const glm::vec3& shadingPoint,
    float u1, float u2, float u3) const {

    if (m_triangles.empty()) {
        return AreaLightSample{};
    }

    // binary search in CDF to select triangle by area
    const auto it = std::lower_bound(m_areaCDF.begin(), m_areaCDF.end(), u1);
    size_t triIndex = std::distance(m_areaCDF.begin(), it);
    triIndex = std::min(triIndex, m_triangles.size() - 1);

    const TriangleData& tri = m_triangles[triIndex];

    // sample point on selected triangle
    AreaLightSample sample = tri.uniformSampler->sample(shadingPoint, u2, u3);

    // adjust pdf for triangle selection
    // pdf_total = pdf_triangle_selection * pdf_point_on_triangle
    // pdf_triangle_selection = area_tri / total_area
    const float selectionProb = tri.area / m_totalArea;
    sample.pdf *= selectionProb;
    sample.area = m_totalArea;
    sample.misWeight = 1.f;

    return sample;
}

float MeshAreaLight::pdfUniform(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    // find which triangle contains the point
    const size_t triIndex = findTriangleContainingPoint(lightPoint);

    if (triIndex >= m_triangles.size()) {
        return 0.0f;  // point not on any triangle
    }

    const TriangleData& tri = m_triangles[triIndex];

    // pdf = (area_tri / total_area) * (1 / area_tri)
    //     = 1 / total_area
    const float trianglePdf = tri.uniformSampler->pdf(shadingPoint, lightPoint);
    const float selectionProb = tri.area / m_totalArea;

    return trianglePdf * selectionProb;
}


/// Area importance sampling
// select triangle proportional to solid angle, then use importance sampling
AreaLightSample MeshAreaLight::sampleAreaImportance(
    const glm::vec3& shadingPoint,
    float u1, float u2, float u3) const {

    if (m_triangles.empty()) {
        return AreaLightSample{};
    }

    // compute solid angle for each triangle
    std::vector<float> solidAngles(m_triangles.size());
    float totalSolidAngle = 0.0f;

    for (size_t i = 0; i < m_triangles.size(); ++i) {
        solidAngles[i] = m_triangles[i].areaImportanceSampler->calculateSolidAngle(shadingPoint);
        if (solidAngles[i] > 0.0f) {
            totalSolidAngle += solidAngles[i];
        }
    }

    // fallback to uniform if no triangle is visible
    if (totalSolidAngle <= 0.0f) {
        return sampleUniform(shadingPoint, u1, u2, u3);
    }

    // build CDF for triangle selection
    std::vector<float> solidAngleCDF(m_triangles.size());
    float cumulative = 0.0f;
    for (size_t i = 0; i < m_triangles.size(); ++i) {
        cumulative += solidAngles[i] / totalSolidAngle;
        solidAngleCDF[i] = cumulative;
    }
    solidAngleCDF.back() = 1.0f;

    // select triangle by solid angle
    const auto it = std::lower_bound(solidAngleCDF.begin(), solidAngleCDF.end(), u1);
    size_t triIndex = std::distance(solidAngleCDF.begin(), it);
    triIndex = std::min(triIndex, m_triangles.size() - 1);

    const TriangleData& tri = m_triangles[triIndex];

    // sample using area importance on selected triangle
    AreaLightSample sample = tri.areaImportanceSampler->sample(shadingPoint, u2, u3);

    // adjust pdf for triangle selection
    const float selectionProb = solidAngles[triIndex] / totalSolidAngle;
    sample.pdf *= selectionProb;
    sample.area = m_totalArea;
    sample.misWeight = 1.f;

    return sample;
}

float MeshAreaLight::pdfAreaImportance(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    // compute solid angles for all triangles
    std::vector<float> solidAngles(m_triangles.size());
    float totalSolidAngle = 0.0f;

    for (size_t i = 0; i < m_triangles.size(); ++i) {
        solidAngles[i] = m_triangles[i].areaImportanceSampler->calculateSolidAngle(shadingPoint);
        if (solidAngles[i] > 0.0f) {
            totalSolidAngle += solidAngles[i];
        }
    }

    if (totalSolidAngle <= 0.0f) {
        return pdfUniform(shadingPoint, lightPoint);
    }

    // find which triangle contains the point
    const size_t triIndex = findTriangleContainingPoint(lightPoint);

    if (triIndex >= m_triangles.size()) {
        return 0.0f;
    }

    const TriangleData& tri = m_triangles[triIndex];

    const float trianglePdf = tri.areaImportanceSampler->pdf(shadingPoint, lightPoint);
    const float selectionProb = solidAngles[triIndex] / totalSolidAngle;

    return trianglePdf * selectionProb;
}


/// Hierarchical flux sampling
// use BVH to select lights proportional to flux/distance²
AreaLightSample MeshAreaLight::sampleHierarchicalFlux(
    const glm::vec3& shadingPoint,
    float u1, float u2, float u3) const {

    // ensure BVH is built
    if (m_bvhNodes.empty() || m_rootNodeIndex == -1) {
        return sampleUniform(shadingPoint, u1, u2, u3);
    }

    // step 1: hierarchically select leaf node
    float pathPdf = 1.0f;
    const int leafNodeIndex = selectBVHLeafNode(shadingPoint, u1, pathPdf);

    if (leafNodeIndex == -1) {
        return AreaLightSample{};
    }

    const BVHNode& leaf = m_bvhNodes[leafNodeIndex];
    if (!leaf.isLeaf || leaf.triangleIndices.empty()) {
        return AreaLightSample{};
    }

    // step 2: select triangle within leaf by flux
    float triangleSelectionProb = 0.0f;
    const size_t selectedTriIdx = selectTriangleFromLeaf(
        leafNodeIndex, u2, triangleSelectionProb);

    if (selectedTriIdx >= m_triangles.size()) {
        return AreaLightSample{};
    }

    const TriangleData& triangle = m_triangles[selectedTriIdx];

    // step 3: sample point on selected triangle
    // reuse random numbers with domain shift to decorrelate
    AreaLightSample sample = triangle.areaImportanceSampler->sample(
        shadingPoint, u3, std::fmod(u1 + u2, 1.0f));

    if (sample.pdf <= 0.0f) {
        return AreaLightSample{};
    }

    // step 4: combine pdfs
    // pdf_total = pdf_path * pdf_triangle_selection * pdf_point
    sample.pdf = pathPdf * triangleSelectionProb * sample.pdf;
    sample.area = m_totalArea;
    sample.misWeight = 1.0f;

    return sample;
}

float MeshAreaLight::pdfHierarchicalFlux(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    if (m_bvhNodes.empty() || m_rootNodeIndex == -1) {
        return pdfUniform(shadingPoint, lightPoint);
    }

    // find which triangle contains the point
    const size_t triIndex = findTriangleContainingPoint(lightPoint);

    if (triIndex >= m_triangles.size()) {
        return 0.0f;
    }

    return computeHierarchicalPdfForTriangle(triIndex, shadingPoint, lightPoint);
}


/// Visibility-aware sampling
// combines hierarchical sampling with learned visibility
AreaLightSample MeshAreaLight::sampleVisibilityAware(
    const glm::vec3& shadingPoint,
    float u1, float u2, float u3) const {

    if (!m_visAwareSampler) {
        // fallback to hierarchical flux if sampler not initialized
        return sampleHierarchicalFlux(shadingPoint, u1, u2, u3);
    }

    // use visibility-aware sampler (includes MIS)
    // note: would need surface normal from hit point in real renderer
    auto sample = m_visAwareSampler->sampleLight(
        shadingPoint,  glm::vec3(0, 1, 0), u1, u2, u3);


    return AreaLightSample{
        sample.position,
        sample.normal,
        sample.pdf,
        sample.misWeight,
        sample.radiance,
        m_totalArea
    };
}


/// BVH management
void MeshAreaLight::buildBVH() {
    if (m_triangles.empty()) {
        return;
    }

    // prepare triangle build info
    std::vector<TriangleBuildInfo> buildInfo;
    buildInfo.reserve(m_triangles.size());

    for (size_t i = 0; i < m_triangles.size(); ++i) {
        const auto& tri = m_triangles[i];

        glm::vec3 bboxMin, bboxMax;
        MathHelpers::computeTriangleBoundingBox(
            tri.v0, tri.v1, tri.v2, bboxMin, bboxMax);

        glm::vec3 centroid = MathHelpers::computeTriangleCentroid(
            tri.v0, tri.v1, tri.v2);

        float flux = MathHelpers::computeTriangleFlux(
            m_emission, tri.intensity, tri.area);

        buildInfo.emplace_back(centroid, bboxMin, bboxMax, flux, tri.area, i);
    }

    // build BVH
    BVHBuilder builder;
    m_rootNodeIndex = builder.build(m_bvhNodes, buildInfo);

    // calculate flux for leaf nodes
    calculateLeafFlux();

    // build mapping from triangles to leaf nodes
    buildTriangleToNodeMapping();

    std::cout << "BVH built with " << m_bvhNodes.size() << " nodes\n";
}

void MeshAreaLight::calculateLeafFlux() {
    for (auto& node : m_bvhNodes) {
        if (node.isLeaf) {
            node.totalFlux = 0.0f;
            for (int triIdx : node.triangleIndices) {
                if (triIdx < static_cast<int>(m_triangles.size())) {
                    const auto& tri = m_triangles[triIdx];
                    node.totalFlux += MathHelpers::computeTriangleFlux(
                        m_emission, tri.intensity, tri.area);
                }
            }
        }
    }
}

void MeshAreaLight::buildTriangleToNodeMapping() {
    m_triangleToLeafNode.resize(m_triangles.size(), -1);

    // traverse BVH and map triangles to their leaf nodes
    std::function<void(int)> traverse = [&](int nodeIndex) {
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_bvhNodes.size())) {
            return;
        }

        const BVHNode& node = m_bvhNodes[nodeIndex];

        if (node.isLeaf) {
            for (int triIdx : node.triangleIndices) {
                if (triIdx < static_cast<int>(m_triangles.size())) {
                    m_triangleToLeafNode[triIdx] = nodeIndex;
                }
            }
        } else {
            traverse(node.leftChild);
            traverse(node.rightChild);
        }
    };

    if (m_rootNodeIndex >= 0) {
        traverse(m_rootNodeIndex);
    }
}

void MeshAreaLight::updateBVH() {
    buildBVH();
}


/// Hierarchical sampling helpers
int MeshAreaLight::selectBVHLeafNode(
    const glm::vec3& shadingPoint,
    float u,
    float& outPathPdf) const {

    outPathPdf = 1.0f;

    if (m_rootNodeIndex < 0 || m_rootNodeIndex >= static_cast<int>(m_bvhNodes.size())) {
        return -1;
    }

    // recursive traversal
    int currentNode = m_rootNodeIndex;

    while (true) {
        const BVHNode& node = m_bvhNodes[currentNode];

        // reached leaf
        if (node.isLeaf) {
            return currentNode;
        }

        // get children
        if (node.leftChild < 0 || node.rightChild < 0 ||
            node.leftChild >= static_cast<int>(m_bvhNodes.size()) ||
            node.rightChild >= static_cast<int>(m_bvhNodes.size())) {
            return -1;
        }

        const BVHNode& left = m_bvhNodes[node.leftChild];
        const BVHNode& right = m_bvhNodes[node.rightChild];

        // calculate importance: flux / distance²
        const glm::vec3 toLeft = left.getCentroid() - shadingPoint;
        float distSqLeft = glm::dot(toLeft, toLeft);
        distSqLeft = std::max(distSqLeft, 1e-4f);

        const glm::vec3 toRight = right.getCentroid() - shadingPoint;
        float distSqRight = glm::dot(toRight, toRight);
        distSqRight = std::max(distSqRight, 1e-4f);

        const float importanceLeft = left.totalFlux / distSqLeft;
        const float importanceRight = right.totalFlux / distSqRight;
        const float totalImportance = importanceLeft + importanceRight;

        if (totalImportance <= 0.0f) {
            return -1;
        }

        // choose child stochastically
        const float probLeft = importanceLeft / totalImportance;

        if (u < probLeft) {
            outPathPdf *= probLeft;
            u = (probLeft > 0.0f) ? u / probLeft : 0.0f;
            currentNode = node.leftChild;
        } else {
            const float probRight = 1.0f - probLeft;
            outPathPdf *= probRight;
            u = (probRight > 0.0f) ? (u - probLeft) / probRight : 0.0f;
            currentNode = node.rightChild;
        }
    }
}

size_t MeshAreaLight::selectTriangleFromLeaf(
    int leafNodeIndex,
    float u,
    float& outSelectionProb) const {

    if (leafNodeIndex < 0 || leafNodeIndex >= static_cast<int>(m_bvhNodes.size())) {
        outSelectionProb = 0.0f;
        return m_triangles.size();
    }

    const BVHNode& leaf = m_bvhNodes[leafNodeIndex];

    if (!leaf.isLeaf || leaf.triangleIndices.empty()) {
        outSelectionProb = 0.0f;
        return m_triangles.size();
    }

    // compute flux for each triangle in leaf
    std::vector<float> fluxes;
    fluxes.reserve(leaf.triangleIndices.size());

    float totalFlux = 0.0f;
    for (int triIdx : leaf.triangleIndices) {
        if (triIdx < static_cast<int>(m_triangles.size())) {
            const auto& tri = m_triangles[triIdx];
            const float flux = MathHelpers::computeTriangleFlux(
                m_emission, tri.intensity, tri.area);
            fluxes.push_back(flux);
            totalFlux += flux;
        }
    }

    if (totalFlux <= 0.0f) {
        outSelectionProb = 0.0f;
        return m_triangles.size();
    }

    // select triangle proportional to flux
    const float target = u * totalFlux;
    float cumulative = 0.0f;
    size_t selectedLocalIdx = 0;

    for (size_t i = 0; i < fluxes.size(); ++i) {
        cumulative += fluxes[i];
        if (target <= cumulative) {
            selectedLocalIdx = i;
            break;
        }
    }

    const int selectedTriIdx = leaf.triangleIndices[selectedLocalIdx];
    outSelectionProb = fluxes[selectedLocalIdx] / totalFlux;

    return static_cast<size_t>(selectedTriIdx);
}

float MeshAreaLight::computeHierarchicalPdfForTriangle(
    size_t triangleIndex,
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    if (triangleIndex >= m_triangleToLeafNode.size()) {
        return 0.0f;
    }

    const int leafNodeIndex = m_triangleToLeafNode[triangleIndex];

    if (leafNodeIndex < 0 || leafNodeIndex >= static_cast<int>(m_bvhNodes.size())) {
        return 0.0f;
    }

    // compute path pdf from root to leaf
    const float pathPdf = computePathPdfToLeaf(leafNodeIndex, shadingPoint);

    if (pathPdf <= 0.0f) {
        return 0.0f;
    }

    // compute triangle selection pdf within leaf
    const BVHNode& leaf = m_bvhNodes[leafNodeIndex];
    float triangleFlux = MathHelpers::computeTriangleFlux(
        m_emission, m_intensity, m_triangles[triangleIndex].area);

    float triangleSelectionProb = 0.0f;
    if (leaf.totalFlux > 0.0f) {
        triangleSelectionProb = triangleFlux / leaf.totalFlux;
    }

    // compute point pdf on triangle
    const float pointPdf = m_triangles[triangleIndex].uniformSampler->pdf(
        shadingPoint, lightPoint);

    return pathPdf * triangleSelectionProb * pointPdf;
}

float MeshAreaLight::computePathPdfToLeaf(
    int leafNodeIndex,
    const glm::vec3& shadingPoint) const {

    // traverse from root to leaf, accumulating probabilities
    float pdf = 1.0f;
    int currentNode = m_rootNodeIndex;

    // helper to check if target leaf is in subtree
    std::function<bool(int, int)> isInSubtree = [&](int target, int root) -> bool {
        if (root == target) return true;
        if (root < 0 || root >= static_cast<int>(m_bvhNodes.size())) return false;

        const BVHNode& node = m_bvhNodes[root];
        if (node.isLeaf) return false;

        return isInSubtree(target, node.leftChild) ||
               isInSubtree(target, node.rightChild);
    };

    while (currentNode != leafNodeIndex) {
        if (currentNode < 0 || currentNode >= static_cast<int>(m_bvhNodes.size())) {
            return 0.0f;
        }

        const BVHNode& node = m_bvhNodes[currentNode];

        if (node.isLeaf) {
            // reached leaf before target - error
            return 0.0f;
        }

        const BVHNode& left = m_bvhNodes[node.leftChild];
        const BVHNode& right = m_bvhNodes[node.rightChild];

        // compute importances (same as in selection)
        const glm::vec3 toLeft = left.getCentroid() - shadingPoint;
        float distSqLeft = std::max(glm::dot(toLeft, toLeft), 1e-4f);

        const glm::vec3 toRight = right.getCentroid() - shadingPoint;
        float distSqRight = std::max(glm::dot(toRight, toRight), 1e-4f);

        const float importanceLeft = left.totalFlux / distSqLeft;
        const float importanceRight = right.totalFlux / distSqRight;
        const float totalImportance = importanceLeft + importanceRight;

        if (totalImportance <= 0.0f) {
            return 0.0f;
        }

        // determine which child leads to target
        const bool goLeft = isInSubtree(leafNodeIndex, node.leftChild);

        if (goLeft) {
            const float prob = importanceLeft / totalImportance;
            pdf *= prob;
            currentNode = node.leftChild;
        } else {
            const float prob = importanceRight / totalImportance;
            pdf *= prob;
            currentNode = node.rightChild;
        }
    }

    return pdf;
}

/// Geometry queries
size_t MeshAreaLight::findTriangleContainingPoint(const glm::vec3& point) const {
    for (size_t i = 0; i < m_triangles.size(); ++i) {
        const auto& tri = m_triangles[i];
        float u, v, w;
        if (isPointInTriangle(point, tri.v0, tri.v1, tri.v2, u, v, w)) {
            return i;
        }
    }
    return m_triangles.size();  // not found
}

bool MeshAreaLight::isPointInTriangle(
    const glm::vec3& point,
    const glm::vec3& v0,
    const glm::vec3& v1,
    const glm::vec3& v2,
    float& u, float& v, float& w) {

    auto bary = MathHelpers::computeBarycentricCoordinates(point, v0, v1, v2);
    u = bary.u;
    v = bary.v;
    w = bary.w;

    return bary.isInside();
}

/// Configuration
void MeshAreaLight::setSamplingStrategy(SamplingStrategy strategy) {
    if (m_strategy == strategy) {
        return;
    }

    m_strategy = strategy;

    // build BVH if needed for hierarchical strategies
    if ((strategy == SamplingStrategy::HierarchicalFlux ||
         strategy == SamplingStrategy::VisibilityAwareHierarchical) &&
        m_bvhNodes.empty()) {
        buildBVH();
    }

    // initialize visibility-aware sampler if needed
    if (strategy == SamplingStrategy::VisibilityAwareHierarchical &&
        s_visibilityCache && !m_visAwareSampler) {

        VisibilityAwareHierarchicalSampler::Configuration config;
        config.enableVisibilityLearning = true;
        config.enableMIS = true;
        config.visibilityWeight = 0.5f;
        config.misHeuristic = MISHeuristic::Balance;

        m_visAwareSampler = std::make_unique<VisibilityAwareHierarchicalSampler>(
            this, s_visibilityCache.get(), config);
    }
}

void MeshAreaLight::setTriangleIntensity(size_t triangleIndex, float intensity) {
    if (triangleIndex >= m_triangles.size()) {
        return;
    }

    TriangleData& tri = m_triangles[triangleIndex];
    tri.intensity = intensity;
    tri.uniformSampler->setIntensity(intensity);
    tri.areaImportanceSampler->setIntensity(intensity);

    // mark BVH as needing rebuild
    // (in production, could do incremental update)
}

/// Accessors
float MeshAreaLight::getTotalPower() const {
    float power = 0.0f;
    for (const auto& tri : m_triangles) {
        power += tri.uniformSampler->getTotalFlux();
    }
    return power;
}

/// Static visibility cache management
void MeshAreaLight::initializeVisibilityCache(
    const glm::vec3& sceneMin,
    const glm::vec3& sceneMax,
    int resolution) {

    s_visibilityCache = std::make_unique<SpatialVisibilityCache>(
        sceneMin, sceneMax, resolution);
}

void MeshAreaLight::getVisibilityStats(
    int& totalCells,
    int& activeCells,
    int& totalSamples) {

    if (s_visibilityCache) {
        s_visibilityCache->getStatistics(totalCells, activeCells, totalSamples);
    } else {
        totalCells = activeCells = totalSamples = 0;
    }
}

void MeshAreaLight::clearVisibilityCache() {
    if (s_visibilityCache) {
        s_visibilityCache->reset();
    }
}

bool MeshAreaLight::containsPoint(const glm::vec3& point, float epsilon) const {
    return findTriangleContainingPoint(point) < m_triangles.size();
}

glm::vec3 MeshAreaLight::getEmissionAt(const glm::vec3& point) const {
    if (containsPoint(point)) {
        return m_emission * m_intensity;
    }
    return glm::vec3(0.0f);
}
