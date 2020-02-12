function [ ] = insertPictureToExcel( filename, picname, sheetnum )
%INSERTPICTURETOEXCEL Summary of this function goes here
%   Detailed explanation goes here

Excel = actxserver('Excel.Application');
Workbooks = Excel.Workbooks;
% Make Excel visible
Excel.Visible=1;
% Open Excel file
Workbook=Workbooks.Open(filename);
% Specify sheet number, data, and range to write to 
Sheets = Excel.ActiveWorkBook.Sheets;
Sheets.Add([], Sheets.Item(Sheets.Count));
sheet = get(Sheets, 'Item', sheetnum);
invoke(sheet, 'Activate');
%Shapes.AddPicture(picname,0,1,400,18,300,235);
sheet.invoke('Pictures').Insert(picname);

invoke(Workbook,'Save')
% Close Excel and clean up
invoke(Excel,'Quit');
delete(Excel);
clear Excel;


end

