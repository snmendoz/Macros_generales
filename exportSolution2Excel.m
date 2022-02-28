function exportSolution2Excel(model,solution,nombreArchivo,met,MetNameFlag)

if nargin<2 || isempty(nombreArchivo); nombreArchivo=['FlujoReacciones_' met];  end
if nargin>=3
    if exist([nombreArchivo '.xls'],'file')
        delete([nombreArchivo '.xls'])
    end
end
if nargin<5; MetNameFlag=1;end
posicionReaccionesActivas=find(solution);
posicionReaccionesActivasSignificativas1=find(solution>1e-6);
posicionReaccionesActivasSignificativas2=find(solution<-1e-6);
posicionReaccionesActivasSignificativas=union(posicionReaccionesActivasSignificativas1,posicionReaccionesActivasSignificativas2);

if MetNameFlag
    %     nombreCompuestos=compartimentalizar(model);
    nombreCompuestos = model.metNames;
else
    nombreCompuestos=model.mets;
end

ecuaciones={};
ecuacionesActivasSignificativas={};
posicionesMetabolitosActivos=[];
posicionesMetabolitosActivosSignificativos=[];

for i=1:length(posicionReaccionesActivas)
    
    %Se encuentran las posiciones de los compuestos de la reaccion i
    PosicionCompuestos=find(model.S(:,posicionReaccionesActivas(i)));
    if isempty(posicionesMetabolitosActivos)
        posicionesMetabolitosActivos=PosicionCompuestos;
    else
        posicionesMetabolitosActivos=union(posicionesMetabolitosActivos,PosicionCompuestos);
    end
    
    
    %Se crea la ecuación que describe la reacción
    reactantes='';
    productos='';
    YaTieneUnReactante=false;
    YaTieneUnProducto=false;
    
    %Se va agregando cada compuesto a la sección de la ecuación que
    %corresponda (reactantes o productos)
    for j=1:length(PosicionCompuestos)
        
        if model.S(PosicionCompuestos(j),posicionReaccionesActivas(i))<0
            if YaTieneUnReactante
                reactantes=[reactantes '  +   ' num2str((model.S(PosicionCompuestos(j),posicionReaccionesActivas(i)))*-1) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
            else
                reactantes=[reactantes num2str((model.S(PosicionCompuestos(j),posicionReaccionesActivas(i)))*-1) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
                YaTieneUnReactante=true;
            end
        else
            if YaTieneUnProducto
                productos=[productos '  +   ' num2str(model.S(PosicionCompuestos(j),posicionReaccionesActivas(i))) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
            else
                productos= [productos num2str(model.S(PosicionCompuestos(j),posicionReaccionesActivas(i))) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
                YaTieneUnProducto=true;
            end
        end
    end
    
    if solution(posicionReaccionesActivas(i))>0
        ecuacion=[reactantes '  ->   ' productos];
    else
        ecuacion=[reactantes '  <-   ' productos];
    end
    
    ecuaciones=[ecuaciones;ecuacion];
end

for i=1:length(posicionReaccionesActivasSignificativas)
    
    %Se encuentran las posiciones de los compuestos de la reaccion i
    PosicionCompuestos=find(model.S(:,posicionReaccionesActivasSignificativas(i)));
    if isempty(posicionesMetabolitosActivos)
        posicionesMetabolitosActivosSignificativos=PosicionCompuestos;
    else
        posicionesMetabolitosActivosSignificativos=union(posicionesMetabolitosActivosSignificativos,PosicionCompuestos);
    end
    
    %Se crea la ecuación que describe la reacción
    reactantes='';
    productos='';
    YaTieneUnReactante=false;
    YaTieneUnProducto=false;
    
    %Se va agregando cada compuesto a la sección de la ecuación que
    %corresponda (reactantes o productos)
    for j=1:length(PosicionCompuestos)
        
        if model.S(PosicionCompuestos(j),posicionReaccionesActivasSignificativas(i))<0
            if YaTieneUnReactante
                reactantes=[reactantes '  +   ' num2str((model.S(PosicionCompuestos(j),posicionReaccionesActivasSignificativas(i)))*-1) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
            else
                reactantes=[reactantes num2str((model.S(PosicionCompuestos(j),posicionReaccionesActivasSignificativas(i)))*-1) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
                YaTieneUnReactante=true;
            end
        else
            if YaTieneUnProducto
                productos=[productos '  +   ' num2str(model.S(PosicionCompuestos(j),posicionReaccionesActivasSignificativas(i))) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
            else
                productos= [productos num2str(model.S(PosicionCompuestos(j),posicionReaccionesActivasSignificativas(i))) ' ' nombreCompuestos{PosicionCompuestos(j)} ' '];
                YaTieneUnProducto=true;
            end
        end
    end
    
    if solution(posicionReaccionesActivasSignificativas(i))>0
        ecuacion=[reactantes '  ->   ' productos];
    else
        ecuacion=[reactantes '  <-   ' productos];
    end
    
    ecuacionesActivasSignificativas=[ecuacionesActivasSignificativas;ecuacion];
end

ecuaciones=['ECUACION';ecuaciones];
flujos=['FLUJO';num2cell(solution(posicionReaccionesActivas))];
Nombres=['NOMBRE REACCION';model.rxnNames(posicionReaccionesActivas)];
xlswrite(nombreArchivo,[ecuaciones flujos Nombres],'R. Activas');

ecuacionesActivasSignificativas=['ECUACION';ecuacionesActivasSignificativas];
flujos=['FLUJO';num2cell(solution(posicionReaccionesActivasSignificativas))];
Nombres=['NOMBRE REACCION';model.rxns(posicionReaccionesActivasSignificativas)];
xlswrite(nombreArchivo,[ecuacionesActivasSignificativas flujos Nombres],'R.A.S.');

IDMetsActivos=['ID METS';model.mets(posicionesMetabolitosActivos)];
NombresMetsActivos=['NAMES METS';model.metNames(posicionesMetabolitosActivos)];
%FormulasMetsActivos=['FORMULAS METS';model.metFormulas(posicionesMetabolitosActivos)];
%xlswrite(nombreArchivo,[IDMetsActivos NombresMetsActivos FormulasMetsActivos],'M. Activas');
xlswrite(nombreArchivo,[IDMetsActivos NombresMetsActivos],'M. Activas');

IDMetsActivos=['ID METS';model.mets(posicionesMetabolitosActivos)];
NombresMetsActivos=['NAMES METS';model.metNames(posicionesMetabolitosActivos)];
%FormulasMetsActivos=['FORMULAS METS';model.metFormulas(posicionesMetabolitosActivos)];
%xlswrite(nombreArchivo,[IDMetsActivos NombresMetsActivos FormulasMetsActivos],'M.A.S.');
xlswrite(nombreArchivo,[IDMetsActivos NombresMetsActivos],'M.A.S.');

disp(['Se imprimieron las reacciones en el archivo: ' nombreArchivo])

end