# Development Manual

## Develop BellComm.jl Pkg

This project is packaged using `Pkg.jl`. For code changes to be reflected in the
packaged software, the package must be set to development mode. This is done by
running the command:

```
julia -e 'using Pkg; Pkg.develop("BellComm")'
```

## Run Tests

Run the BellComm.jl package tests (continuous integration):

```
julia -e 'using Pkg;  Pkg.test("BellComm")'
```

Download test dependencies into local environment:

```
julia --project=test/ -e 'using Pkg; Pkg.instantiate()'
```

Run tests from a dev environment:

```
julia test/path/to/test_file.jl
```

## Build Docs

Deploy docs server on local machine:

```
cd docs/build/; python3 -m http.server --bind localhost
```

Build docs from BellComm.jl source (continuous integration):

1. Download docs dependencies into local environment:

        julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'


2. Build docs:

        julia --project=docs/ docs/make.jl

Build docs from a dev environment:

```
julia -e 'using Pkg; Pkg.add("Documenter")'
```

```
julia docs/make.jl
```
