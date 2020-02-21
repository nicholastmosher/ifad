# `ifad`

Interactive Functional Annotations Dashboard

## Installation

In order to run this project, you'll need to install Rust. Head on over to
[rustup.rs] to find installation instructions. On
Linux/Mac, this should be a one-line terminal command. On Windows, you'll
just download a `.exe` file, but you'll also need to make sure you have the
[Microsoft Visual C++ compiler] installed. After you install Rust, you may
need to restart your terminal for the new commands to show up.

[rustup.rs]: https://rustup.rs/
[Microsoft Visual C++ compiler]: https://visualstudio.microsoft.com/vs/features/cplusplus/

## Getting Started

This repository contains a command-line tool to quickly play with some of the
Gene and Annotation filtering options. We can use `cargo run` to quickly see
what options are available:

```
$ cargo run --release
ifad
USAGE:
    ifad [OPTIONS]
FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
OPTIONS:
        --genes=<genes>                        the file to read genes from (e.g. gene-types.txt
        --annotations=<annotations>            the file to read annotations from (e.g. tair.gaf)
        --genes-out=<genes_out>                the file to write queried genes to (e.g. gene-types_F-EXP.txt
        --annotations-out=<annotations_out>    the file to write queried annotations to (e.g. tair_F-EXP.gaf)
        --query=<query>                        the type of query [default: union]  [possible values: union,
                                               intersection]
        --segment=<segment>...                 a segment to use in the query, given as ASPECT,STATUS (e.g. F,EXP or
                                               C,OTHER)%
```

We can see that there are a few required options. The first ones to notice are
`--genes` and `--annotations`, which let you choose the source data files to
query over. Next are the `--genes-out` and `--annotations-out` options, which
let you choose filenames to save the queried subsets of data into. Finally,
there's the `--segment` option, which lets you choose "segments" of the data
you'd like to have exported. You can use the `--segment` option multiple times.

Let's try it out. I'll assume you have a `gene-types.txt` file and a `tair.gaf`
file ready.

```
$ cargo run --release -- \
            --genes=./gene-types.txt \
            --annotations=./tair.gaf \
            --genes-out=./gene-types_F-EXP.txt \
            --annotations-out=./tair_F-EXP.gaf \
            --segments=F,EXP
```

After running this, you'll notice two new files have been created,
`gene-types_F-EXP.txt` and `tair_F-EXP.gaf`, with the subsets of gene data
and annotation data that belong to `F,EXP` (Molecular Function with Experimental Evidence).
