#!/bin/sh

# make skymap movies
ffmpeg -y -i Frames/map_ModelA_%d_sm20.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p ModelA.mp4
ffmpeg -y -i Frames/map_ModelB_%d_sm20.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p ModelB.mp4
ffmpeg -y -i Frames/map_ModelC_%d_sm20.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p ModelC.mp4
ffmpeg -y -i Frames/map_ModelB_no4945_%d_sm20.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p ModelB_no4945.mp4

# make composition movies
ffmpeg -y -i Frames/lightmap_ModelC_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p light_ModelC.mp4
ffmpeg -y -i Frames/heavymap_ModelC_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p heavy_ModelC.mp4
ffmpeg -y -i Frames/lightmap_ModelB_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p light_ModelB.mp4
ffmpeg -y -i Frames/heavymap_ModelB_%d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p heavy_ModelB.mp4

# make timeseries movies
ffmpeg -y -i Frames/summary_ModelA_%03d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p summary_ModelA.mp4
ffmpeg -y -i Frames/summary_ModelB_%03d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p summary_ModelB.mp4
ffmpeg -y -i Frames/summary_ModelC_%03d.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p summary_ModelC.mp4

# make spatial distribution movies
#ffmpeg -y -i Frames/xy_ModelA_%03dd.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p xy_ModelA.mp4
#ffmpeg -y -i Frames/xy_ModelB_%03dd.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p xy_ModelB.mp4
#ffmpeg -y -i Frames/xy_ModelC_%03dd.png -filter:v "setpts=3*PTS" -pix_fmt yuv420p xy_ModelC.mp4
