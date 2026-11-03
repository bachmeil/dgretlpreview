# Build and run test programs on Windows (LDC, MSVC target).
# Usage: testing\build_win.ps1 [test_name ...]     (default: every test_*.d here)
#
# Needs: ldc2 and Visual Studio Build Tools (the linker); libxml2, zlib, OpenBLAS and FFTW3 installed with
#   vcpkg install libxml2 zlib openblas clapack fftw3 --triplet x64-windows   (the vcpkg OpenBLAS has no LAPACK; clapack supplies it)
# Set VCPKG_ROOT, or have vcpkg in %USERPROFILE%\repos\vcpkg.
#
# Each source file is compiled to its own object file (in parallel): compiling all of them in one ldc2 run crashed LDC 1.41 on
# Windows (inside std.format; perhaps only because vcpkg's bin was on PATH, see below). Objects are rebuilt only when the source is newer.
param([string[]]$Tests, [switch]$NoRun)

$ErrorActionPreference = 'Stop'
$here = $PSScriptRoot
$root = Split-Path $here
$out = Join-Path $here 'build'
$obj = Join-Path $out 'obj'
New-Item -ItemType Directory -Force $obj | Out-Null

$vcpkg = if ($env:VCPKG_ROOT) { $env:VCPKG_ROOT } else { Join-Path $env:USERPROFILE 'repos\vcpkg' }
$vlib = Join-Path $vcpkg 'installed\x64-windows\lib'
$vbin = Join-Path $vcpkg 'installed\x64-windows\bin'
# The vcpkg DLLs (zlib, iconv, libxml2) go on PATH only while a test runs: with them on PATH ldc2 itself loads them
# instead of its own copies and crashes (hence the "random" compiler crashes seen with LDC 1.41).
$origPath = $env:PATH
# test_gnuplot_spawn runs gnuplot when GNUPLOT_EXE names it (it is skipped otherwise)
if (-not $env:GNUPLOT_EXE) {
    $gp = (Get-Command gnuplot -ErrorAction SilentlyContinue).Source
    if (-not $gp -and (Test-Path "$env:ProgramFiles\gnuplot\bin\gnuplot.exe")) { $gp = "$env:ProgramFiles\gnuplot\bin\gnuplot.exe" }
    if ($gp) { $env:GNUPLOT_EXE = $gp }
}

function Build-Objects($sources) {
    $todo = $sources | Where-Object {
        $o = Join-Path $obj (($_.FullName.Substring($root.Length + 1) -replace '[\\/]', '_') -replace '\.d$', '.obj')
        -not (Test-Path $o) -or (Get-Item $o).LastWriteTime -lt $_.LastWriteTime
    }
    $todo | ForEach-Object -Parallel {
        $o = Join-Path $using:obj (($_.FullName.Substring($using:root.Length + 1) -replace '[\\/]', '_') -replace '\.d$', '.obj')
        $msg = ldc2 -c -g -mscrtlib=msvcrt -I"$using:root" "-of=$o" $_.FullName 2>&1 | Out-String
        if ($LASTEXITCODE -ne 0) { Write-Host "FAILED $($_.Name)`n$msg" }
    } -ThrottleLimit 8
    $sources | ForEach-Object { Join-Path $obj (($_.FullName.Substring($root.Length + 1) -replace '[\\/]', '_') -replace '\.d$', '.obj') }
}

# The dgretl modules are compiled together: their generated functions (postblits of private structs, for example) are
# only emitted into the object file of a module that is compiled in the same run as the one that uses them.
# (No -g here: debug info makes LDC 1.41 crash on this group.)
function Build-Group($sources) {
    $stamp = Join-Path $obj 'dgretl.stamp'
    $newest = ($sources | Measure-Object LastWriteTime -Maximum).Maximum
    $objs = $sources | ForEach-Object { Join-Path $obj (($_.BaseName -eq 'package') ? 'dgretl.obj' : ('dgretl.' + $_.BaseName + '.obj')) }
    if (-not (Test-Path $stamp) -or (Get-Item $stamp).LastWriteTime -lt $newest) {
        $msg = ldc2 -c -oq -mscrtlib=msvcrt -I"$root" "-od=$obj" ($sources | ForEach-Object FullName) 2>&1 | Out-String
        if ($LASTEXITCODE -ne 0) { Write-Host "FAILED dgretl`n$msg" } else { Set-Content $stamp 'ok' }
    }
    $objs
}
# The import library that vcpkg makes for openblas.dll also defines printf, sprintf, ... (the DLL exports its own copies
# of the CRT's inline stdio functions), which hijack the C runtime's. Make a clean one with only the BLAS entry points.
function New-CleanOpenBlas {
    $dll = Join-Path $vbin 'openblas.dll'
    $clean = Join-Path $out 'openblas_clean.lib'
    if ((Test-Path $clean) -and (Get-Item $clean).LastWriteTime -gt (Get-Item $dll).LastWriteTime) { return $clean }
    $vsroot = & "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -property installationPath
    $msvc = Get-ChildItem "$vsroot\VC\Tools\MSVC" | Sort-Object Name | Select-Object -Last 1
    $bin = Join-Path $msvc.FullName 'bin\HostX64\x64'
    $names = & "$bin\dumpbin.exe" /exports $dll | ForEach-Object {
        if ($_ -match '^\s+\d+\s+[0-9A-F]+\s+[0-9A-F]{8}\s+(\S+)') { $Matches[1] } } |
        Where-Object { $_ -match '_$|^cblas_|^openblas_' }
    $def = Join-Path $out 'openblas_clean.def'
    @('LIBRARY openblas', 'EXPORTS') + $names | Set-Content $def
    $env:PATH = "$bin;$env:PATH"    # lib.exe needs the linker's neighbours
    & "$bin\lib.exe" /nologo /machine:x64 "/def:$def" "/out:$clean" | Out-Null
    $clean
}
# Run ldc2 with a time limit; sets $LASTEXITCODE. LDC 1.41 on Windows crashes now and then ("illegal instruction")
# on the same command that works the next time, so a failed run is repeated (up to $tries times).
function Invoke-Ldc($ldcArgs, $seconds = 300, $tries = 3) {
    foreach ($i in 1..$tries) {
        $p = Start-Process ldc2 -ArgumentList ($ldcArgs | ForEach-Object { '"' + $_ + '"' }) -NoNewWindow -PassThru
        if (-not $p.WaitForExit($seconds * 1000)) {
            $p.Kill($true)
            Write-Host "ldc2 timed out after $seconds s"
            $global:LASTEXITCODE = 1
        } else {
            $global:LASTEXITCODE = $p.ExitCode
        }
        if ($global:LASTEXITCODE -eq 0) { return }
    }
}
$skip = 'slowPrng.d', 'testgslprng.d'
$lib = Get-ChildItem "$root\gretl\*.d", "$root\gretl\plugin\*.d", "$root\gretl\minpack\*.d", "$root\cephes\*.d" |
    Where-Object { $_.Name -notin $skip }
$dgretl = Get-ChildItem "$root\dgretl\*.d" | Where-Object { $_.Name -notin $skip }
$stubs = Get-Item "$here\stubs.d"
$unported = Get-Item "$here\stubs_unported.d"

if (-not $Tests) { $Tests = Get-ChildItem "$here\test_*.d" | ForEach-Object { $_.BaseName } }

$cleanBlas = New-CleanOpenBlas
$libObjs = Build-Objects ($lib + $stubs + $unported)
$dgObjs = $null
$status = 0
$failed = @()
foreach ($t in $Tests) {
    Write-Host "=== building $t"
    $src = Get-Item "$here\$t.d"
    $objs = $libObjs
    if (Select-String -Path $src -Pattern '^import dgretl' -Quiet) {
        if (-not $dgObjs) { $dgObjs = Build-Group $dgretl }
        $objs = $libObjs + $dgObjs
    }
    $exe = Join-Path $out "$t.exe"
    $tobj = Join-Path $obj "test_$t.obj"
    Invoke-Ldc @('-c', '-mscrtlib=msvcrt', "-I$root", "-of=$tobj", $src.FullName)
    if ($LASTEXITCODE -ne 0) {
        # LDC 1.41 on Windows crashes while compiling some of the larger tests (inside std.format's floating-point code);
        # without asserts and contracts (-release) they compile. Only the test program is affected, not the library.
        Write-Host "retrying $t with -release"
        Invoke-Ldc @('-c', '-release', '-mscrtlib=msvcrt', "-I$root", "-of=$tobj", $src.FullName)
    }
    if ($LASTEXITCODE -ne 0) { $status = 1; $failed += "$t (build)"; continue }
    # (absolute paths, and objects only: LDC 1.41 crashes on some other combinations)
    # The vcpkg libraries are DLL builds, so the C runtime is the dynamic one. Phobos already contains zlib.
    ldc2 -mscrtlib=msvcrt "-L/NODEFAULTLIB:libcmt" "-of=$exe" $tobj $objs `
        (Join-Path $vlib libxml2.lib) $cleanBlas (Join-Path $vlib lapack.lib) (Join-Path $vlib libf2c.lib) (Join-Path $vlib fftw3.lib) (Join-Path $vlib iconv.lib) (Join-Path $vlib charset.lib)
    if ($LASTEXITCODE -ne 0) { $status = 1; $failed += "$t (build)"; continue }
    if (-not $NoRun) {
        Write-Host "=== running $t"
        Push-Location $here    # the tests use paths relative to testing\
        $env:PATH = "$vbin;$origPath"
        & $exe
        $code = $LASTEXITCODE
        $env:PATH = $origPath
        Pop-Location
        if ($code -ne 0) { $status = 1; $failed += "$t (exit $code)" }
    }
}
Write-Host "=== $($Tests.Count) tests, $($failed.Count) failed"
$failed | ForEach-Object { Write-Host "  $_" }
exit $status
