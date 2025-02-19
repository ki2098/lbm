import struct
import numpy as np

with open('TurbSim.bts', mode='rb') as f:
    id, nz, ny, ntwr, nt = struct.unpack('<h4l', f.read(2+4*4))
    dz, dy, dt, uhub, zhub, zbottom = struct.unpack('<6f', f.read(6*4))
    scl = np.zeros(3, np.float32)
    off = np.zeros(3, np.float32)
    scl[0], off[0], scl[1], off[1], scl[2], off[2] = struct.unpack('<6f', f.read(6*4))
    nchar, = struct.unpack('<l', f.read(4))
    info = f.read(nchar).decode()

    print(id, nz, ny, ntwr, nt)
    print(dz, dy, dt, uhub, zhub, zbottom)
    print(scl, off)

    u = np.zeros((3, nt, ny, nz))
    utwr = np.zeros((3, nt, ntwr))

    for it in range(nt):
        buffer = np.frombuffer(
            f.read(2*3*ny*nz), dtype=np.int16
        ).astype(np.float32).reshape([3, ny, nz], order='F')
        u[:, it, :, :] = buffer

        buffer = np.frombuffer(
            f.read(2*3*ntwr), dtype=np.int16
        ).astype(np.float32).reshape([3, ntwr], order='F')
        utwr[:, it, :] = buffer
    
    print(u[0,1,1,1])

    u -= off[:, None, None, None]
    u /= scl[:, None, None, None]
    utwr -= off[:, None, None]
    utwr /= scl[:, None, None]

    print(u[0,1,1,1])

    for it in range(0, nt, 20):
        fname = 'data/uvw.csv.%d'%(it)
        with open(fname, "w") as fo:
            fo.write('x,y,z,u,v,w\n')
            for k in range(nz):
                for j in range(ny):
                    fo.write(
                        '%f,%f,%f,%f,%f,%f\n'%(
                            0, dy*j, dz*k, u[0,it,j,k], u[1,it,j,k], u[2,it,j,k]
                        )
                    )