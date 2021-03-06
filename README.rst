Spectral Archetypes
===================

A partial repository of code related to spectral modeling.  Methods and code may include the following:

* Heteroskedastic Matrix Factorization (Tsalmantza & Hogg)
* Integer-programming-generated archetypes (Roweis & Hogg)
* Platykurtic vector discovery (Hogg)

Authors
-------

* **David W. Hogg** (NYU)
* **Vivi Tsalmantza** (MPIA)
* **Sam Roweis** (deceased)
* **Jo Bovy** (IAS)

License
-------

Copyright 2009, 2010, 2011, 2012, 2013, 2014 the authors.  **All rights reserved**

Notes
-----

Migrated from SVN at <http://Astrometry.net/> by the following::

    cd
    git svn clone svn+ssh://astrometry.net/svn/trunk/projects/archetypes/ \
        --no-minimize-url --authors-file ~/authors
    cd archetypes
    git pull git@github.com:davidwhogg/SpectralArchetypes.git
    git remote add origin git@github.com:davidwhogg/SpectralArchetypes.git
    git push origin master

(and something similear for `svn+ssh://astrometry.net/svn/trunk/documents/papers/archetypes/`.)
