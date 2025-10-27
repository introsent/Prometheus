//
// Created by minaj on 10/4/2025.
//

#include "embree_scene.h"
EmbreeScene::EmbreeScene(const EmbreeDevice* devicePtr) {
    m_scene = rtcNewScene(devicePtr->handle());
    rtcSetSceneBuildQuality(m_scene, RTC_BUILD_QUALITY_HIGH);
}

void EmbreeScene::commit() {
    rtcCommitScene(m_scene);
}

void EmbreeScene::release() {
    rtcReleaseScene(m_scene);
}

RTCScene EmbreeScene::handle() const {
    return m_scene;
}