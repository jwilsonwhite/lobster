function Connmat = load_connmat(type)

% type 1 = average
% type 3 = individual years

%if type == 1
%Connmat = importdata('connectivity matrix years1+2+3 lobster oct28 2016-v2.txt');
%Connmat = Connmat'; % transpose to obtain sources on columns, destinations on rows
%elseif type == 2
for i = 1:3
  %  x = importdata(strcat('year',num2str(i),'-lobster.txt'));
    x = importdata(strcat('year',num2str(i),'-lobster apr24 no grounds.txt'));
    Connmat(:,:,i) = x';
end
%end

if type == 1
    Connmat = mean(Connmat,3);
end
    