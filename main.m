     
%reset variables
clear variables

%read setup

for i = 3:3

    filename = ['setup\setup',num2str(i),'.txt'];

    setfile = importdata(filename);
    for setnumber = 1:size(setfile.data,1)
        theset = setfile.data(setnumber,:);
        testpop(setnumber,theset,i);
    end
end

% setfile2 = importdata('setup2.txt');
% for setnumber = 1:size(setfile2,1)
%     E_testpop(setnumber,setfile2(setnumber,:));
% end