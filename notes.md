# Para Cargar ROOT y coatjava

Estos son los pasos para cargar ROOT y de coatjava

## 1) OPTIONAL (Solo si alguna vez te dice "module: Command not found")
    ```
    source /etc/profile.d/modules.csh
    ```    

## 2) Decirle dónde están los módulos de JLab en CVMFS
    ```
    module use /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/modulefiles/ module load clas12
    ```

## 3) OPTIONAL Ver qué ROOT hay disponible
    ```
    module avail root
    ```

## 4) ROOT ya se carga por defecto desde el module clas12
    Para verificar cual módulo de ROOT está cargado
    ```
    root-config --version
    which root
    ```

## 5) Ver qué coatjava hay disponible
    ```
    module avail coatjava
    ```

## 6) Cargar el modulo de coatjava que deseamos
    ```
    module load coatjava/13.5.2
    ```

## 7) Verificar cual hipo-utils se está usando
    ```
    which hipo-utils
    ```

## 8) Ver un solo archivo hipo
    ```
    hipo-utils -dump /work/clas12b/users/lixu/singleParticle/testOSG/output0/dst-0.hipo
    ```

## 9) Quiero guardar la información de los hipo files de las tres ramas
    Debo abrir la terminal y luego poner ```bash``` para entrar en un entorno de bash
    ```
    bash

    out="/w/hallb-scshelf2102/clas12/nlbucuru/CNDHitFinder_Optimization/hipo_dump_all.txt"
    : > "$out"

    base="/work/clas12b/users/lixu/singleParticle"
    folders=("testOSG" "testCJ0" "testCJ1")

    for name in "${folders[@]}"; do
    for i in {0..19}; do
        f="$base/$name/output$i/dst-$i.hipo"
        if [[ -f "$f" ]]; then
        echo "==================== $f ====================" >> "$out"
        hipo-utils -info "$f" >> "$out"
        echo >> "$out"
        fi
    done
    done

    echo "Wrote: $out"

    ```

.L analysis_DVCS_preskimmed_fiducials_NID.C+
.L spring2019_preskimmed.C+
.L run_DVCS.C+


run_DVCS("pDVCS", "spring2019");  
run_DVCS("pDVCS", "fall2019");   
run_DVCS("pDVCS", "spring2020"); 
