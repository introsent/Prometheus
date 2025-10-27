//
// Created by minaj on 10/4/2025.
//

#include "embree_device.h"

EmbreeDevice::EmbreeDevice() {
    m_device = rtcNewDevice(NULL);
}

void EmbreeDevice::release() {
    rtcReleaseDevice(m_device);
}

RTCDevice EmbreeDevice::handle() const {
    return m_device;
}


