#!/bin/bash
start=$(date +%s)
export DOCKER_BUILDKIT=1
# 使用 Docker buildx 的计时功能
docker buildx build \
    --progress=plain \
    --build-arg BUILDKIT_INLINE_CACHE=1 \
    --network=host \
    -t deephall_v4  \
    -f ~/workspace/code/DeepHall/DeepHallZeros/docker/build_docker/Dockerfile_zj_v4.deephall \
    ~/workspace/code/DeepHall

#    --cache-from type=local,src=/tmp/buildx-cache \
#     --cache-to type=local,dest=/tmp/buildx-cache-new \


end=$(date +%s)
runtime=$((end-start))

echo "Build took: $runtime seconds"
echo "Build took: $((runtime/60)) minutes and $((runtime%60)) seconds"

#     -t deephall_v1 \
#     --target netobs_installer \