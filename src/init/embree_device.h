//
// Created by minaj on 10/4/2025.
//

#ifndef EMBREE_INIT_H
#define EMBREE_INIT_H
#include <embree4/rtcore.h>
class EmbreeDevice {
public:
    EmbreeDevice();

    void release();
    [[nodiscard]] RTCDevice handle() const;
private:
    RTCDevice m_device = nullptr;
};
#endif //EMBREE_INIT_H
