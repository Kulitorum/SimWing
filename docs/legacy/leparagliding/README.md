# Inherited LEparagliding engineering records

These documents were active planning, release, and subsystem records in the
LEparagliding Studio repository at the SimWing bootstrap commit recorded in
`UPSTREAM.md`. They are archived here to keep SimWing's active documentation
focused on the coupled XPBD-CFD program.

They are not disposable history. The imported editor, exact-model, Print/Cut,
Playground, and soft-wing code still embodies decisions and invariants recorded
here. Read the applicable record before modifying those inherited subsystems:

- `CONTINUE.md`: current inherited Playground/free-flight handoff and measured
  regression guards;
- `playground-shape-analysis.md`: shape/metric rationale and claim boundary;
- `xpbd-performance.md`: structural solver and aerodynamic load-path history;
- `BACKLOG.md`: deliberately deferred editor and Print/Cut work;
- `flat-part-orientation.md`: unresolved fabric-grain engineering problem;
- `CLAUDE.md`: historical machine, release, manual, and signing procedures.

The active SimWing direction is
[`../../coupled-fsi-architecture.md`](../../coupled-fsi-architecture.md). When
an inherited subsystem is replaced, migrate any surviving constraint into an
active SimWing design record before removing its archived source.
