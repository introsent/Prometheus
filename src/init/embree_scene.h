//
// Created by minaj on 10/4/2025.
//
#ifndef EMBREE_SCENE_H
#define EMBREE_SCENE_H

#include "../init/embree_device.h"
class EmbreeScene {
public:
    explicit EmbreeScene(const EmbreeDevice* devicePtr);

    void commit();
    void release();
    [[nodiscard]] RTCScene handle() const;

private:
    RTCScene m_scene = nullptr;
};
#endif //EMBREE_SCENE_H
