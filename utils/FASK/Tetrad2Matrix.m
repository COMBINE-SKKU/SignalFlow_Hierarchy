function [mat]=Tetrad2Matrix(tet,d) %,outputfile)
% TETRAD2MATRIX transform a list of Tetrad output edges into a matrix form
%
%INPUT
%*tet -- tetrad ouput .txt file containing list of nodes followed by list of edges
%*d   -- if d='undirected' make a binary UNDIRECTED SYMMETRIC MATRIX, 
%        if d='directed' make a binary DIRECTED MATRIX. 
%           directions goes row to column i --> j (for Smith's simulation)
%        if d='pattern' encode 2 for --> and 1 for ---       
%       
%
%*outputfile -- name for the output .txt file containing the matrix
%
%OUPUT
% A txt file containing the matrix, in tab delimited format
%

%tic

%read the file with the list of edges and the name of the variables.
%tet=importdata(tet); one function but very very slow. 

%read the file and save each line in a cell row, exactly as importdata but
%way way faster
fid=fopen(tet);
tline = fgetl(fid); %each execution of fgetl reads the next line
tet = cell(0,1);  %define the cell to save each line
while ischar(tline)  %while the line was characters
    tet{end+1,1} = tline;  %iteratively save the lines in the cel
    tline = fgetl(fid); %run fgetl to read the next line
end
%tet(end) = []; %Tetrad saves an extra line with no characters that messes things up. Remove it here.
fclose(fid);


%create a list of all the nodes. This is done to assign the exact position of the edges in the matrix.
w{1}='Graph Nodes:';
a=find(strcmp(tet,w))+1; %find the beginning of the nodes list
v{1}='Graph Edges:';  %if not working add a space after the colon :
b=find(strcmp(tet,v))-1; %find the end of the nodes list
D=[tet{a:b}]; %put together all the names of the nodes
R=textscan(D,'%s','Delimiter',';'); %split them into separate elements if csv
%R=textscan(D,'%s'); %if labels separated by space
R=R{:}; %transform them into a list of strings

%determine where the list of edges start in the tetrad txt file
q=find(strcmp(tet,v))+1;

e{1}='Ambiguous'; 
f=find(strncmp(tet,e,3))-2;

if isempty(f)
    fi=length(tet);
else
    fi=f;
end

%make the matrix nxn where n is the number of nodes.
mat=zeros(length(R),length(R));


for i=q:fi               %length(tet)
  
     n=strsplit(tet{i}); 
      
     %do not consider the entries of motion variables
         if strncmp(n{2},'Motion',6) || strncmp(n{4},'Motion',6)
             %jump to next case
         
         else
     
     %node 1
  %      [~,endIndex] = regexp(n{2},'VAR'); 
  %      n1 = str2double(n{2}(endIndex+1:end)); 
     %node 2
  %      [~,endIndex] = regexp(n{4},'VAR'); 
  %      n2 = str2double(n{4}(endIndex+1:end));
   
     
%      %if the labels are not X# then use this
%      %node 1
      n1=find(strcmp(R,n{2})); % assign the position in the matrix of the corrsponding AAL roi
%      %node 2
      n2=find(strcmp(R,n{4})); % assign the position in the matrix of the corrsponding AAL roi

     %arrowtail
     at=sscanf(n{3},'%c%*c%*c');
     %arrowhead
     ah=sscanf(n{3},'%*c%*c%c');
     

     
    if strcmp(d,'undirected'); 
        %make an undirected symmetric matrix.
        mat(n1,n2)=1;
        mat(n2,n1)=1;
    
    elseif strcmp(d,'directed');
        %make a binary directed matrix, 
        if strcmp(at,'-') && strcmp(ah,'>')   %|| strcmp(at,'-') && strcmp(ah,'-')
           mat(n1,n2)=1;
        end
        
    elseif strcmp(d,'pattern');
        %encode directed edges as 2, and undirected edges as 1
        if strcmp(at,'-') && strcmp(ah,'>')
            mat(n1,n2)=2;
        elseif strcmp(at,'-') && strcmp(ah,'-')
            mat(n1,n2)=1;
        end
        
        
    end
         end
         
end