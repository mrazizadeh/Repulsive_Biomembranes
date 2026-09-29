@echo off
setlocal
REM Launch Repulsive Biomembranes via WSL2 + WSLg without attaching a console.
REM Attaching to cmd.exe lets Windows "copy mode" freeze the GUI.
REM Optional argument: scene file inside test_case (default: sceneADE.txt)

set SCENE=%~1
if "%SCENE%"=="" set SCENE=sceneADE.txt

REM Wake WSL / WSLg so Weston can attach a real monitor before the app starts.
wsl.exe -d Ubuntu-24.04 -u khaled -- true >nul 2>&1

powershell.exe -NoProfile -WindowStyle Hidden -Command ^
  "Start-Process -FilePath 'wsl.exe' -WindowStyle Hidden -ArgumentList '-d','Ubuntu-24.04','-u','khaled','--','bash','-lc','set +u; source /opt/intel/oneapi/setvars.sh --force >/dev/null; export DISPLAY=:0 WAYLAND_DISPLAY=wayland-0 XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir; cd /home/khaled/Repulsive_Biomembranes/test_case; exec /home/khaled/Repulsive_Biomembranes/build/bin/biorsurfaces \"%SCENE%\" > /tmp/biorsurfaces.log 2>&1'"

echo Launched Repulsive BioSurfaces. Look for that title on the taskbar.
echo Log: \\wsl$\Ubuntu-24.04\tmp\biorsurfaces.log
endlocal
