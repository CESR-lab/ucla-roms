# Obtaining and compiling ROMS

## Obtaining ROMS
Begin by cloning this repository to a suitable location with

```git clone https://github.com/CWorthy-ocean/ucla-roms.git```

You will need to set the ``ROMS_ROOT`` environment variable to the location of your clone:

```bash
cd ucla-roms
export ROMS_ROOT=$(pwd)
```

As with any environment variables described on the [previous page](dependencies), `ROMS_ROOT` needs to be set every time you wish to use ROMS, and it may be convenient to set its value in a location such as a shell startup file (e.g. `~/.bashrc`)

## Compiling internal dependencies

The ROMS repository includes two additional directories that should be compiled the first time ROMS is installed:

```bash
cd $ROMS_ROOT/NHMG/src
make
```

compiles the optional non-hydrostatic modeling library, while

```bash
cd $ROMS_ROOT/Tools-Roms/
make
```

compiles a series of programs used to manage ROMS input and output. More on these later.

