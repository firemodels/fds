% McDermott
% 9-28-2011
% svn_log_report.m

close all
clear all
fclose all

% read xml file created from (for example):
%
% rmcdermo@blaze FDS-SMV]$ svn log --xml >> fds-smv_log.xml

tagName = 'author'
findLabel = 'randy.mcdermott'
to_revision = 6951;

xDoc = xmlread('fds-smv_log.xml')

% xmlwrite(xDoc)
%
%   <logentry revision="7">
%
%      <author>bryanwklein</author>
%
%      <date>2007-05-01T21:17:11.111165Z</date>
%
%      <msg>Created directory FDS/trunk</msg>
%
%   </logentry>

allListitems = xDoc.getElementsByTagName('logentry');

fid = fopen('log_output.txt','wt');

for k = 0:allListitems.getLength-1
    
    thisListitem = allListitems.item(k);
    
    thisList = thisListitem.getElementsByTagName(tagName);
    thisElement = thisList.item(0);
    
    if strcmp(thisElement.getFirstChild.getData, findLabel)
        
        revision = thisListitem.getAttribute('revision');
        if str2num(char(revision))<to_revision; break; end
        display(['revision = ',char(revision)]);
        fprintf(fid,'%s\n',['revision = ',char(revision)]);
        
        thisList = thisListitem.getElementsByTagName('author');
        thisElement = thisList.item(0);
        author = char(thisElement.getFirstChild.getData);
        display(['author = ',author]);
        fprintf(fid,'%s\n',['author = ',author]);
        
        thisList = thisListitem.getElementsByTagName('date');
        thisElement = thisList.item(0);
        date = char(thisElement.getFirstChild.getData);
        display(['date = ',date]);
        fprintf(fid,'%s\n',['date = ',date]);
        
        thisList = thisListitem.getElementsByTagName('msg');
        thisElement = thisList.item(0);
        msg = char(thisElement.getFirstChild.getData);
        display(['msg = ',msg]);
        fprintf(fid,'%s\n',['msg = ',msg]);
        
        display([' ']); % blank line
        fprintf(fid,'%s\n',[' ']);
    end
    
end

fclose(fid);

