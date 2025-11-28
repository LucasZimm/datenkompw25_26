@echo off
setlocal enabledelayedexpansion

set PROJECT_ROOT=%~dp0
set BIN_DIR=%PROJECT_ROOT%bin\MSVC-19.42.34436.0\Release
set INPUT_DIR=%PROJECT_ROOT%kodak-pgm
set OUTPUT_ROOT=%PROJECT_ROOT%output

set CSV_FILE=%OUTPUT_ROOT%\bpp_results.csv

if not exist "%OUTPUT_ROOT%" mkdir "%OUTPUT_ROOT%"

echo Datei;Codec;Breite;Hoehe;Pixel;Bin_Bytes;BPP > "%CSV_FILE%"

REM Alle EXEs durchlaufen
for %%E in (%BIN_DIR%\*.exe) do (
    set "EXE_PATH=%%~fE"
    set "EXE_NAME=%%~nE"
    set "EXE_OUTPUT_DIR=%OUTPUT_ROOT%\!EXE_NAME!"

    if not exist "!EXE_OUTPUT_DIR!" mkdir "!EXE_OUTPUT_DIR!"

    echo ▶ Starte: !EXE_NAME!

    if not exist "!EXE_PATH!" (
        echo   EXE nicht gefunden: !EXE_PATH!
        goto :continueExeLoop
    )

    REM Alle .pgm Dateien verarbeiten
    for %%F in (%INPUT_DIR%\*.pgm) do (
        set "INPUT_FILE=%%~fF"
        set "ENCODE_FILE=%%~nF.bin"
        set "DECODE_FILE=%%~nF.pgm"

        echo   Kodieren: %%~nxF
        "!EXE_PATH!" -e "!INPUT_FILE!" "!ENCODE_FILE!" "!EXE_OUTPUT_DIR!"

        echo   Dekodieren: !ENCODE_FILE!
        "!EXE_PATH!" -d "!EXE_OUTPUT_DIR!\!ENCODE_FILE!" "!DECODE_FILE!" "!EXE_OUTPUT_DIR!"

        REM Dateigröße der BIN holen
        for %%S in ("!EXE_OUTPUT_DIR!\!ENCODE_FILE!") do set FILE_SIZE=%%~zS
        echo   Dateigröße Bytes !FILE_SIZE!

        REM Breite und Höhe aus PGM lesen (2. Zeile mit Zahlen)
        set WIDTH=
        set HEIGHT=
        set linenum=0
        for /f "usebackq tokens=1,2" %%a in ("!INPUT_FILE!") do (
            set /a linenum+=1
            if !linenum! equ 2 (
                set WIDTH=%%a
                set HEIGHT=%%b
            )
        )
        echo Breite=!WIDTH!, Höhe=!HEIGHT!

        REM Pixel berechnen
        set /a PIXELS=WIDTH*HEIGHT
        set /a BITS=FILE_SIZE*8
        echo   Breite !WIDTH!, Höhe !HEIGHT!, Pixel !PIXELS!

        REM BPP mit 3 Nachkommastellen berechnen via PowerShell
        for /f %%B in ('powershell -NoProfile -Command "Write-Output ([math]::Round((!BITS! / !PIXELS!),3))"') do set BPP=%%B
        echo   BPP: !BPP!

        REM Ergebnis ins CSV schreiben
        echo %%~nxF;!EXE_NAME!;!WIDTH!;!HEIGHT!;!PIXELS!;!FILE_SIZE!;!BPP! >> "%CSV_FILE%"
    )
    :continueExeLoop
    echo.
)

echo ✅ Fertig! Ergebnis: %CSV_FILE%
pause
