import numpy as np
# import pandas as pd


# A simulation data

t1=np.arange(1, 50, step=1)
t2=np.arange(1, 49, step=1)

u1=np.sin(t1*np.pi*2/49)
v1=np.cos(t2*np.pi*2/48)

signal=np.outer(u1, v1)


data=signal+np.random.default_rng().normal(0, 0.2, size=(49, 48))




X = np.arange(1, 49, 1)
Y = np.arange(1, 50, 1)
X, Y = np.meshgrid(X, Y)



def svd3dplot(data, ncomp=3, isurface=True, *args):
    import numpy as np
    from matplotlib import pyplot as plt
    nrow, ncol=data.shape
    u, s, vh=np.linalg.svd(data)
    v=vh.transpose()
    umat=u[:, :ncomp]
    svec=s[:ncomp]
    vmat=v[:, :ncomp]
    app=np.matmul(umat, np.diag(svec))
    app=np.matmul(app, vmat.transpose())
    res=data-app
    compmat=np.empty([nrow, ncol, ncomp])
    for i in np.arange(ncomp):
        utemp=u[:, i]
        stemp=s[i]
        vtemp=v[:, i]
        compmat[:, :, i]=stemp*np.outer(utemp, vtemp)
        
    #Now we will start to draw plots
    X = np.arange(1, ncol+1, 1)
    Y = np.arange(1, nrow+1, 1)
    X, Y = np.meshgrid(X, Y)   
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    
    if isurface:
        fig = plt.figure(figsize=(8, 12))
    #    ax = fig.add_subplot(2, 3, 1, projection='3d')
        ax=plt.subplot(321, projection='3d')
        surf = ax.plot_surface(X, Y, data, cmap="jet",
                               linewidth=1, antialiased=True)

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.3, aspect=5)
        ax.title.set_text('Original Data')
        
    #    ax = fig.add_subplot(2, 3, 2, projection='3d')
        ax=plt.subplot(322, projection='3d')

        surf = ax.plot_surface(X, Y, compmat[:, :, 0], cmap="jet",
                               linewidth=1, antialiased=True)

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.3, aspect=5)
        ax.title.set_text('Component 1')
        
    #    ax = fig.add_subplot(2, 3, 3, projection='3d')
        ax=plt.subplot(323, projection='3d')
        surf = ax.plot_surface(X, Y, compmat[:, :, 1], cmap="jet",
                               linewidth=1, antialiased=True)

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.3, aspect=5)
        ax.title.set_text('Component 2')

        ax=plt.subplot(324, projection='3d')
    #    ax = fig.add_subplot(2, 3, 4, projection='3d')
        surf = ax.plot_surface(X, Y, compmat[:, :, 2], cmap="jet",
                               linewidth=1, antialiased=True)

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.3, aspect=5)
        ax.title.set_text('Component 3')
        
        ax=plt.subplot(325, projection='3d')
    #    ax = fig.add_subplot(2, 3, 2, projection='3d')
        surf = ax.plot_surface(X, Y, app, cmap="jet",
                               linewidth=1, antialiased=True)

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.3, aspect=5)
        ax.title.set_text('Approximation')

        ax=plt.subplot(326, projection='3d')
    #   ax = fig.add_subplot(2, 3, 2, projection='3d')
        surf = ax.plot_surface(X, Y, res, cmap="jet",
                               linewidth=1, antialiased=True)

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.3, aspect=5)
        ax.title.set_text('Residual')
        plt.show()

    if (not isurface):
        fig = plt.figure(figsize=(8, 12))

        ax=plt.subplot(321)
        im=plt.imshow(data, cmap="jet")
        plt.colorbar(im, shrink=0.85)
        ax.title.set_text('Original Data')

        ax=plt.subplot(322)
        im=plt.imshow(compmat[:, :, 0], cmap="jet")
        plt.colorbar(im, shrink=0.85)
        ax.title.set_text('Component 1')
    
        ax=plt.subplot(323)
        im=plt.imshow(compmat[:, :, 1], cmap="jet")
        plt.colorbar(im, shrink=0.85)
        ax.title.set_text('Component 2')

        ax=plt.subplot(324)
        im=plt.imshow(compmat[:, :, 2], cmap="jet")
        plt.colorbar(im, shrink=0.85)
        ax.title.set_text('Component 3')
    
        ax=plt.subplot(325)
        im=plt.imshow(app, cmap="jet")
        plt.colorbar(im, shrink=0.85)
        ax.title.set_text('Approximation')
    
        ax=plt.subplot(326)
        im=plt.imshow(res, cmap="jet")
        plt.colorbar(im, shrink=0.85)
        ax.title.set_text('Residual')
        plt.show()



#### showing the examples

## Image view with 3 components
svd3dplot(data, isurface=False)

## Surface view with 3 components
svd3dplot(data)
    