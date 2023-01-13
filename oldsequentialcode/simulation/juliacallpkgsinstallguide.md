# To install packages which can be used on juliacall do the following:

1) Go to your environment's julia folder and open manifest.toml
    - search for package's name, if you find it then copy its uuid 
    --> if the package is not there, install it: 
        - open Julia REPL
        - go into pkg (press ])
        - run (juliaenv) pkg> add <pkgname>
        - let installation go
        - now the uuid will be found in manifest.toml

2) Open a python script either python or jupyter 
    - import juliapkg
    - run the code: juliapkg.add(<pkgname>, uuid = <uuid from manifest.toml>)

3) Now you can import the package into your juliacall code like:
    import juliacall as jl
    jl.seval("using packagename")

4) GG 