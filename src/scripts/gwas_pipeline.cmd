REM Script loops through .xml configFiles for TASSEL.

for %%f in (rs_*.xml) do (
    echo %%~nf
    run_pipeline.bat -configFile "%%f"
    )