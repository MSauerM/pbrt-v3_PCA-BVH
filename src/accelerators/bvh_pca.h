//
// Created by msauer on 12.06.19.
//

#ifndef PBRT_V3_BVH_PCA_H
#define PBRT_V3_BVH_PCA_H

#include <core/primitive.h>

namespace pbrt{
    struct BVHBuildNode;

    struct BVHPrimitiveInfo;

    struct PCAPrimitiveInfo;

    struct LinearBVHNode;


    class PCAAccel: public Aggregate{
        public:
        enum class SplitMethod {SAH, Middle};

        PCAAccel(std::vector<std::shared_ptr<Primitive>> p, int maxPrimsInNode = 1,
                SplitMethod splitMethod = SplitMethod::Middle);

        Bounds3f WorldBound() const;
        ~PCAAccel();
        bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
        bool IntersectP(const Ray &ray) const;

        private:
        BVHBuildNode *recursiveBuild(MemoryArena &arena,
                std::vector<BVHPrimitiveInfo> &primitiveInfo,
                int start, int end, int *totalNodes,
                std::vector<std::shared_ptr<Primitive>> &orderedPrims);

        int flattenBVHTree(BVHBuildNode *node, int *offset);

        const int maxPrimsInNode;
        const SplitMethod splitMethod;
        std::vector<std::shared_ptr<Primitive>> primitives;
        LinearBVHNode *nodes = nullptr;
    };


    std::shared_ptr<PCAAccel> CreatePCAAccelerator(
            std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

    float CalculateCost(BVHBuildNode* currentNode);



}


#endif //PBRT_V3_BVH_PCA_H
