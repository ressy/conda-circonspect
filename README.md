# conda-circonspect

A simple wrapper around [Circonspect] to install it to a dedicated [Anaconda]
environment.  (Ideally this would be an Anaconda package itself, but for now,
just a wrapper.)  The bundled Circonspect code is the latest stable release
listed on Sourceforge as of this writing.

    ./setup.sh
    source activate circonspect
    Circonspect -h

[Circonspect]: https://sourceforge.net/projects/circonspect/
[Anaconda]: https://anaconda.org/
