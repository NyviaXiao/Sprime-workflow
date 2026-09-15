# Compile the bundled map_arch only when config points to the bundled binary.
# An externally managed executable is treated as an ordinary input instead.
if Path(C['map_arch']).resolve() == (Path(ROOT) / 'tools/map_arch/map_arch').resolve():
    rule compile_map_arch:
        input: source=str(Path(ROOT) / 'tools/map_arch/map_arch.v2.c'), makefile=str(Path(ROOT) / 'tools/map_arch/makefile')
        output: C['map_arch']
        params: stage='compile', code=CODE
        log: 'logs/compile_map_arch.log'
        script: TASK
