using MovingFrame
using GeometryTypes
using CodecZlib
using FileIO
fname = "/home/yusri/pCloudDrive/190612/0429/gc-50_0.2_1_0.2/each/t0001/cell00001.tif.gz"
shape = GZip.open(fname) do io
    FileIO.load(io)
end
