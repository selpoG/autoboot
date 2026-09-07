# GAPToMathematica

## Setup

Install the [.NET 10 SDK](https://dotnet.microsoft.com/download/dotnet/10.0)
and [GAP](https://www.gap-system.org/). GNU Make 4.3 or later is also required.
On macOS, install GNU Make with Homebrew and replace `make` with `gmake` in
the commands below.

Then run:

```sh
make
```

This will generate `GAPToMathematica.dll` in `bin/`.

## Usage

```sh
dotnet bin/GAPToMathematica.dll
```

Alternatively, build and run it in one step:

```sh
make run
```
