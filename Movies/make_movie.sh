#!/bin/sh

ffmpeg -y -i Frames/map_ModelB_80-20_%d_sm20.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p ModelB.mp4
ffmpeg -y -i Frames/map_ModelC_80-20_%d_sm20.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p ModelC.mp4
ffmpeg -y -i Frames/lightmap_ModelC_80-20_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p light_ModelC.mp4
ffmpeg -y -i Frames/heavymap_ModelC_80-20_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p heavy_ModelC.mp4
ffmpeg -y -i Frames/lightmap_ModelB_80-20_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p light_ModelB.mp4
ffmpeg -y -i Frames/heavymap_ModelB_80-20_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p heavy_ModelB.mp4