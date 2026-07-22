# Commit Convention

Just to keep everything a bit more clean.

Use Conventional Commits:

```text
<type>(<scope>): <summary>
```

Examples:

```text
feat(dwt): implement forward Haar transform
fix(core): prevent out-of-bounds convolution
test(operators): add odd-length coverage
docs(readme): update project status
build(makefile): add sanitizer target
```

## Types

| Type | Use |
|---|---|
| `feat` | Add user-visible or public API functionality |
| `fix` | Correct a defect |
| `refactor` | Restructure code without changing behavior |
| `perf` | Improve performance |
| `test` | Add or update tests |
| `docs` | Change documentation only |
| `build` | Change the build system or dependencies |
| `ci` | Change continuous integration |
| `chore` | Perform repository maintenance |
| `revert` | Revert an earlier commit |
