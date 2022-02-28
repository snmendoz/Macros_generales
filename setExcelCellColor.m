function setExcelCellColor(fileName,cell, sheet,color)

Excel = actxserver('excel.application');
% Get Workbook object
WB = Excel.Workbooks.Open(fullfile(pwd, fileName),0,false);
% Set the color of cell "A1" of Sheet 1 to RED (3)
WB.Worksheets.Item(sheet).Range(cell).Interior.ColorIndex = color;
% Save Workbook
WB.Save();
% Close Workbook
WB.Close();
% Quit Excel
Excel.Quit();

end