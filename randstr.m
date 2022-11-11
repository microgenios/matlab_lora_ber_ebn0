function strchar = randstr(dim,dim3,varargin)
%% Purpose:
% ----------
%  This routine will generate a random string sequence consisting of the
%  same number of characters and/or numbers as specified by the dim input
%
%  Example function call:
%
%  exStr = randstr([2 5],2,'useWildCards',false);
%
%% Inputs:
%  ----------
%
%  dim                     [N x M x P x  ...]       Character array size
%                                                   specifying the
%                                                   dimensions of the
%                                                   return character string
%                                                   sequence
%
%  dim3                     integer                 Specifier of singleton
%                                                   dimension in which the
%                                                   string is formed
%                                                   e.g. when dim = [1 10]
%                                                   the singleton dimension
%                                                   is dim3 = 2 such that
%                                                   the string reads across 
%                                                   columns. 
%
%
% Optional Input Arguments:
% ---------------------------
%
%  useDigits                boolean/integer         if boolean,
%                                                   {true} = use numbers in
%                                                   the randomly generated
%                                                   string
%
%                                                   false = do not use
%                                                   numbers in randomly 
%                                                   generated 
%                                                   string
%
%                                                   integer specifies how
%                                                   many numbers to use in
%                                                   the randomly generated
%                                                   string of length
%                                                   dim(dim3)
%
% useWildCards              boolean/integer         if boolean,
%                                                   {true} = use wild card 
%                                                   characters in the  
%                                                   randomly generated
%                                                   string
%
%                                                   false = do not use
%                                                   wild card characters in 
%                                                   the randomly generated 
%                                                   string
%
%                                                   integer specifies how
%                                                   many wild card characters 
%                                                   to use in the randomly
%                                                   generated string of 
%                                                   length dim(dim3)
%
%
% useLowerCase              boolean/integer         if boolean,
%                                                   {true} = use lower case
%                                                   characters in the  
%                                                   randomly generated
%                                                   string
%
%                                                   false = do not use
%                                                   lower case characters in 
%                                                   the randomly generated 
%                                                   string
%
%                                                   integer specifies how
%                                                   many lower case characters 
%                                                   to use in the randomly
%                                                   generated string of 
%                                                   length dim(dim3)
%
%
% useUpperCase              boolean/integer         if boolean,
%                                                   {true} = use upper case
%                                                   characters in the  
%                                                   randomly generated
%                                                   string
%
%                                                   false = do not use
%                                                   upper case characters in 
%                                                   the randomly generated 
%                                                   string
%
%                                                   integer specifies how
%                                                   many upper case characters 
%                                                   to use in the randomly
%                                                   generated string of 
%                                                   length dim(dim3)
%
%% Outputs:
% ----------
%
%  strchar                  [N x M x P x  ...]      Randomly generated
%                                                   chracter string of
%                                                   dimensions
%                                                   corresponding to dim
%
%% Revision History:
% -----------------
%  Darin C. Koblick                                         (c) 02/24/2020
%% ----------------------  Begin Code Sequence ----------------------------
if nargin == 0
    dim = [2 6 3];
   dim3 = 2;
   strchar = randstr(dim,dim3,'useDigits',3,'useWildCards',false,'useUpperCase',false);
   return;
end
      strchar = NaN;
    useDigits = true;
 useWildCards = true;
 useLowerCase = true;
 useUpperCase = true;
     varNames = {'useDigits','useWildCards','useLowerCase','useUpperCase'};
%% Parse the varaiable input arguments for specific properties:
idx = find(strcmpi(varargin,varNames{1}));
if ~isempty(idx)
       useDigits = varargin{idx+1};
end
idx = find(strcmpi(varargin,varNames{2}));
if ~isempty(idx)
    useWildCards = varargin{idx+1};
end
idx = find(strcmpi(varargin,varNames{3}));
if ~isempty(idx)
    useLowerCase = varargin{idx+1};
end
idx = find(strcmpi(varargin,varNames{4}));
if ~isempty(idx)
    useUpperCase = varargin{idx+1};
end
%% Argument Name Error Check:
          idxStr = cellfun(@isstr,varargin);
          idxVar = ismember(varargin(idxStr),varNames);
if sum(idxVar) ~= sum(idxStr)
        allVars = varargin(idxStr);
    unknownVars = allVars(~idxVar);
    fprintf(2,'%s',[mfilename, ...
             ':: Error: Uknown Variable(s): ',unknownVars{1}]);
    for tuv=2:numel(unknownVars)
       fprintf(2,'%s',[',',unknownVars{tuv}]);
    end
    fprintf(2,'\n'); 
     return;
end          
%% Input Formatting:    
         totNumChar = prod(dim);
          totNumStr = totNumChar/dim(dim3);
totNumCharPerString = totNumChar/totNumStr;
          charCount = 0;                    
%% Error Checking:           
if isnumeric(useDigits)
    charCount = charCount+useDigits;
end
if isnumeric(useWildCards)
    charCount = charCount+useWildCards;
end
if isnumeric(useLowerCase)
    charCount = charCount+useLowerCase;
end
if isnumeric(useUpperCase)
    charCount = charCount+useUpperCase;
end
%% Quantity Error Check
if (charCount < dim(dim3) && charCount > 0 && (~islogical(useDigits) && ...
   ~islogical(useLowerCase) && ~islogical(useUpperCase) &&  ...
   ~islogical(useWildCards))) || charCount > dim(dim3)
   fprintf(2,'%s\n',[mfilename, ...
             ':: Error: Sum of Specified Quantities Must Equal ', ...
             num2str(dim(dim3))]);
   return;
end
%% Initialize sub-matricies:
      randDigits = [];
   randLowerCase = [];
   randUpperCase = [];
    randWildCase = [];
       randASCII = [];       
%% Create full ASCII lookup table:
charCodes = 33:126;
    ascii = cellfun(@(x)char(x),num2cell(charCodes),'UniformOutput',false);
%% Determine subcases of ascii:
   asciiDigits = ascii(~cellfun(@isempty,regexp(ascii,'[0-9]')));
asciiLowerCase = ascii(~cellfun(@isempty,regexp(ascii,'[a-z]')));    
asciiUpperCase = ascii(~cellfun(@isempty,regexp(ascii,'[A-Z]'))); 
asciiWildCards = ascii(~ismember(ascii,[asciiDigits,asciiLowerCase,asciiUpperCase]));
%Formulate the string character options for the unspecified amounts:
   asciiTable = [];
   if useDigits && islogical(useDigits)
       asciiTable = [asciiTable, asciiDigits];
   end
   if useLowerCase && islogical(useLowerCase)
       asciiTable = [asciiTable, asciiLowerCase];
   end
   if useUpperCase && islogical(useUpperCase)
       asciiTable = [asciiTable, asciiUpperCase]; 
   end
   if useWildCards && islogical(useWildCards)
       asciiTable = [asciiTable, asciiWildCards]; 
   end
%% With the remaining count, draw from the ascii table:
if ~isempty(asciiTable)
    totNumASCII = (totNumCharPerString-charCount)*totNumStr;
      randASCII = reshape(asciiTable(randi(numel(asciiTable), ...
                          totNumASCII,1)),[totNumStr (totNumCharPerString-charCount)]);  
end   
%% Based on the quantities specified by the user, perform a random draw:    
   if ~islogical(useDigits) && useDigits > 0
            totNumDigits = useDigits*totNumStr;
              randDigits = reshape(asciiDigits(randi(numel(asciiDigits), ...
                                   totNumDigits,1)),[totNumStr useDigits]);     
   end
   if ~islogical(useLowerCase) && useLowerCase > 0
         totNumLowerCase = useLowerCase*totNumStr;
           randLowerCase = reshape(asciiLowerCase(randi(numel(asciiLowerCase), ...
                                   totNumLowerCase,1)),[totNumStr useLowerCase]);
   end
   if ~islogical(useUpperCase) && useUpperCase > 0
         totNumUpperCase = useUpperCase*totNumStr;
           randUpperCase = reshape(asciiUpperCase(randi(numel(asciiUpperCase), ...
                                   totNumUpperCase,1)),[totNumStr useUpperCase]);
   end
   if ~islogical(useWildCards) && useWildCards > 0
         totNumWildCards = useWildCards*totNumStr;
            randWildCase = reshape(asciiWildCards(randi(numel(asciiWildCards), ...
                                   totNumWildCards,1)),[totNumStr useWildCards]);
   end
%% Concatinate the array and shuffle columns:
    strchar = [randASCII,randDigits,randLowerCase,randUpperCase,randWildCase];
   for tc=1:totNumStr
        strchar(tc,:) = strchar(tc,randperm(totNumCharPerString));
   end
%% Convert to a matrix:
        strchar = cell2mat(strchar);
%% Reshape the string array such that it matches the original dimension:
            [~,fSeq] = fDim(NaN(dim),dim3);
             strchar = eDim(strchar,fSeq);     
end

function eND = eDim(fND,fSeq)
%% Purpose:
%Take any flattened 2-D matrix in MATLAB and convert it back into its
%original multi-dimensional form before flatening. This is often times 
%necessary when writing complex operations on multi-dimensional matrices.  
%It is also desired that after flattening, the dimension that is preserved 
%has the correct sequence. This is especially important for vector 
%processing. Once flattened, and an operation has been performed on the 2-D
%matrix, often times the 2-D matrix will need to be converted back to the
%original multi-dimensional matrix.  This function will extract the
%multi-dimensional matrix from the flattened 2-D matrix using fDim().
%
%% Inputs:
% fND           [N x size(ND,dim)]                  A 2-D matrix with
%                                                   N corresponding to the
%                                                   same number of elements
%                                                   in the ND matrix but
%                                                   with a single dimension
%                                                   preserved.
%
% fSeq           struct                             fSeq.dim is the initial
%                                                   location of the desired
%                                                   dimension(s) before the
%                                                   shift sequence.
%                                                   
%                                                   fSeq.postShift is the
%                                                   matrix dimensions
%                                                   before the reshaping
%                                                   sequence.
%
%
%% Outputs:
% eND            [O x P x Q x R x S x ...]          Unflattened matrix with 
%                                                   any number and order of
%                                                   dimensions.
%
%
%% Created By Darin Koblick (C) 07/19/2012
iND = reshape(fND,fSeq.postShift);
%inject a singleton dimension:
if numel(size(iND)) ~= numel(fSeq.postShift)
    %shift dimensions to the right and add singletons along the way:
    eND = shiftdim(iND,numel(size(iND))-numel(fSeq.postShift));
    %Add a circular shift to the left and wrap the leading edges:
    eND = shiftdim(eND,numel(fSeq.postShift)-(fSeq.dim)+1);
else
    %Shift dimensions to the left and wrap the leading edges:
    eND = shiftdim(iND,numel(fSeq.postShift)-(fSeq.dim));
end
end

function [fND,fSeq] = fDim(ND,dim)
%% Purpose:
%Take any N-D matrix in MATLAB and flatten it down into an N x size(ND,dim)
%2-D matrix.  This is often times necessary when writing complex operations 
%on multi-dimensional matrices.  It is also desired that after flattening, 
%the dimension that is preserved has the correct sequence. This is especially
%important for vector processing.
%
%% Inputs:
% ND            [O x P x Q x R x S x ...]           Any matrix with any
%                                                   number, N, of
%                                                   dimensions.
%
% dim           double                              Specify a single
%                                                   dimension in which to
%                                                   preserve when creating
%                                                   the columns of the 2-D
%                                                   flattened matrix.
%
%% Outputs:
% fND           [N x size(ND,dim)]                  A 2-D matrix with
%                                                   N corresponding to the
%                                                   same number of elements
%                                                   in the ND matrix but
%                                                   with a single dimension
%                                                   preserved.
%
% fSeq           struct                             fSeq.dim is the initial
%                                                   location of the desired
%                                                   dimension(s) before the
%                                                   shift sequence.
%                                                   
%                                                   fSeq.postShift is the
%                                                   matrix dimensions
%                                                   before the reshaping
%                                                   sequence.
%% Created By Darin Koblick (C) 07/19/2012
%Take the desired dimension and shift it up-front before reshaping the array
if dim > ndims(ND)
    dim = 0;
end 
iND = shiftdim(ND,dim);
fSeq.postShift = size(iND);
fSeq.dim = dim;
if dim > 0
    fND = reshape(iND,numel(ND)/size(ND,dim),size(ND,dim));
else
    fND = reshape(iND,numel(ND),1);
end
end