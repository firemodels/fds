#!/usr/bin/perl
use LWP::Simple;
#
$yearl="y:";
$dayl="d:";
$hourl="h";
$mailto="glenn.forney\@nist.gov kevin.mcgrattan\@nist.gov william.mell\@nist.gov";
@email_alert=("glenn.forney\@nist.gov", "gforney1\@verizon.net", "kevin.mcgrattan\@nist.gov", "bryan.klein\@nist.gov");
$cell_alert="3018070458\@vtext.com";
$cell2_alert="3012528200\@messaging.nextel.com";
#$cell_alert="glenn.forney\@nist.gov";
@hosts=("fire70","fire71","fire72","fire73","fire74","fire75","fire76","fire77","fire78","fire79","fire61","fire63","fire64","fire65","fire66","fire67","fire68","fire62","fire51","fire52","fire53","fire54","fire55","fire56","fire57","fire58","fire59","fire42","fire43","fire44","fire45","fire46","fire47","leo","fire32","fire33","fire34","fire35","capella","hadar");
%hostnames=("capella","fire11/capella","arcturus","fire12/arcturus","spica","fire13/spica","hadar","fire14/hadar","leo","fire31/leo","fire32","fire32","fire33","fire33","fire34","fire34","fire35","fire35","acrux","fire41/acrux","fire42","fire42","fire43","fire43","fire44","fire44","fire45","fire45","fire46","fire46","fire47","fire47","fire51","fire51","fire52","fire52","fire53","fire53","fire54","fire54","fire55","fire55","fire56","fire56","fire57","fire57","fire58","fire58","fire59","fire59","fire61","fire61","fire62","fire62","fire63","fire63","fire64","fire64","fire65","fire65","fire66","fire66","fire67","fire67","fire68","fire68","fire70","fire70","fire71","fire71","fire72","fire72","fire73","fire73","fire74","fire74","fire75","fire75","fire76","fire76","fire77","fire77","fire78","fire78","fire79","fire79",);
%OS=("capella","RH73","arcturus","RH73","speca","RH73","hadar","RH73","leo","F6","fire32","F6","fire33","F6","fire34","F6","fire35","F6","acrux","F6","fire42","F6","fire43","F6","fire44","F6","fire45","F6","fire46","F6","fire47","F6","fire51","F6","fire52","F6","fire53","F6","fire54","F6","fire55","F6","fire56","F6","fire57","F6","fire58","F6","fire59","F6","fire61","F6","fire62","F6","fire63","F6","fire64","F6","fire65","F6","fire66","F6","fire67","F6","fire68","F6","fire70","F8","fire71","F8","fire72","F8","fire73","F8","fire74","F8","fire75","F8","fire76","F8","fire77","F8","fire78","F8","fire79","F8");
%hostmem=("capella",2048,"spica",2048,"arcturus",2048,"hadar",2048,"leo",2048,"fire32",2048,"fire33",2048,"fire34",2048,"fire35",2048,"acrux",2048,"fire42",4096,"fire43",4096,"fire44",2048,"fire45",2048,"fire46",2048,"fire47",2048,"fire51",2048,"fire52",2048,"fire53",2048,"fire54",2048,"fire55",2048,"fire56",2048,"fire57",2048,"fire58",2048,"fire59",2048,"fire61",,4096,"fire62",8192,"fire63",4096,"fire64",4096,"fire65",4096,"fire66",4096,"fire67",4096,"fire68",4096,"fire70",8192,"fire71",8192,"fire72",8192,"fire73",8192,"fire74",8192,"fire75",8192,"fire76",8192,"fire77",8192,"fire78",8192,"fire79",8192);
%cpuspeed=("leo","3.2","fire32","3.2","fire33","3.2","fire34","3.2","fire35","3.2","capella",2.4,"arcturus","2.4","spica","2.4","hadar","2.4","acrux","3.6","fire42","3.6","fire43","3.6","fire44","3.6","fire45",3.6,"fire46",3.6,"fire47",3.6,"fire51",3.8,"fire52",3.8,"fire53",3.8,"fire54",3.8,"fire55",3.8,"fire56",3.8,"fire57",3.8,"fire58",3.8,"fire59",3.8,"fire61",5.2,"fire62",4.8,"fire63",5.2,"fire64",5.2,"fire65",5.2,"fire66",5.2,"fire67",5.2,"fire68",5.2,"fire70",6.5,"fire71",6.5,"fire72",6.5,"fire73",6.5,"fire74",6.5,"fire75",6.5,"fire76",6.5,"fire77",6.5,"fire78",6.5,"fire79",6.5);
@month_days=(31,28,31,30,31,30,31,31,30,31,30,31);

$statdir="/home/fire51/gforney/STATS/";
$warnfile=$statdir . "WARN_sent_warning_email_high_temp";
$temperature_file=$statdir . "temp_history.csv";
$cellfonefile=$statdir . "WARN_called_cell_fone";
$passwdfile=$statdir . "passwd.yp";
$logdir=$statdir . "LOGS/";
$lockfile=$logdir . "genweb_lock";
$notifylockfile=$logdir . "notify_genweb_lock";
#
# determine if another copy of this script is running now
#
if(-e $lockfile){
# another copy is running, so abort after seeing if we need
# to send a mail message warning of a problem
#
  if(! -e $notifylockfile){
    open(DATE,"date|");
    $lockdate=<DATE>;
    close(DATE);
    open(MAIL,"|mail -s \"SYSTEM STATUS: genwepage.pl script did not complete\" $mailto");
    print MAIL <<EOF;
The script genwebpage.pl did not complete at $lockdate.  A system may be down,
causing this script to hang.
EOF
    close(MAIL);
    system("touch $notifylockfile");
  }
  exit;
}
system("touch $lockfile");
#
# determine which systems are up
#
$nnodes=0;
foreach $host (@hosts){
  $statusfile = $logdir . $host . "_status"; 
  open(STATUSFILE,"<$statusfile");
  $line=<STATUSFILE>;
  close(STATUSFILE);
  chomp($line);
  $systemstatus{"$host"}=$line;
  if($systemstatus{"$host"} != 0){next;}
  $nnodes=$nnodes+2;
}
$sf="<font size=-1>";
$sfc="</font>";
#
# get cluster room temperature
#
$clustertemp="Unavailable";
$doc = get 'http://129.6.160.15/temp';
@temps=split('\|',$doc);
$degsymbol = ' <sup>o</sup>';

$clustertempval=@temps[1];

$degF = $clustertempval;
$degF = int(4*$degF+0.5)/4.0;
$degF = $degF. "<sup>o</sup>" ."F";

$degC = $clustertempval-32;
$degC = 5*$degC/9;
$degC = int($degC+0.5);
$degC = $degC. "<sup>o</sup>" ."C";

#
# cell phone called (once) when temperature exceeds 100F
#
if($clustertempval>95&& !-e $cellfonefile){
  system("touch $cellfonefile");
  open(MAIL,"|mail -s \"*** Clust Temp=$degF, too high\" $cell_alert");
  print MAIL <<EOF;
cluster room, 224/B355, temp=$degF too high.
Call glenn forney @ 301 807 0458 to open door.
EOF
close(MAIL);
  open(MAIL,"|mail -s \"*** Clust Temp=$degF, too high\" $cell2_alert");
  print MAIL <<EOF;
cluster room, 224/B355, temp=$degF too high.
Call glenn forney @ 301 807 0458 to open door.
EOF
close(MAIL);
}
#
# emails sent when temperature exceeds 90 F
#
if($clustertempval>90&& !-e $warnfile){
  system("touch $warnfile");
  foreach $email (@email_alert){
    open(MAIL,"|mail -s \"*** Clust Temp=$degF, too high\" $email");
  print MAIL <<EOF;
***Cluster room (224/B339) temperature=($degF) is too high.
Please open door ASAP.

Questions call:
Glenn Forney
301 975 2313 (w)
301 972 2233 (h)
301 807 0458 (c)

Kevin McGrattan
301 975 2712 (w)

Bryan Klein
301 975 5171 (w)
EOF
close(MAIL);
  }
}

#
# get date, times
#
open(DATETIME,"date +\"%l:%M %p %Z\"|");
$ctime=<DATETIME>;
chomp($ctime);
close(DATETIME);
#
open(DATE,"date +%j|");
$date=<DATE>;
chomp($date);
close(DATE);
#
open(DATE,"date +\"%B %d, %Y\"|");
$date2=<DATE>;
chomp($date2);
close(DATE);
#
open(DATE,"date +\"%B %d\"|");
$date3=<DATE>;
chomp($date3);
close(DATE);
#
open(DATE,"date +\"%B %Y\"|");
$txtmonth=<DATE>;
chomp($txtmonth);
close(DATE);
#
open(DATE,"date +%m|");
$month=<DATE>;
chomp($month);
close(DATE);
#
open(DATE,"date +%k|");
$hour=<DATE>;
chomp($hour);
close(DATE);
#
open(DATE,"date +%M|");
$minute=<DATE>;
chomp($minute);
close(DATE);
 
open(DATE,"date +%j|");
$dayofyear=<DATE>;
chomp($dayofyear);
close(DATE);
#
#
open(DATE,"date +%d|");
$day=<DATE>;
chomp($day);
close(DATE);
$dec_hour = ((1+$minute)/60.0+$hour)/24.0;
$dec_date = ($day-1+((1+$minute)/60.0+$hour)/24.0)/30.0;

$dec_date2 = $dayofyear-1+($minute/60.0+$hour)/24.0;
system("echo $dec_date2 $clustertempval>>$temperature_file");
system("cp $statdir/temp_history.png /var/www/html/load/.");

system($statdir . "makesummaries.pl");
#
open(PASSWD,"ypcat passwd|");
open(PASSWDOUT,">$passwdfile");
$line=<PASSWD>;
while($line ne ""){
  chomp($line);
  print PASSWDOUT ("$line\n");
  ($user,$dummy,$uid)=split(/:/,$line);
  $users{"$uid"}=$user;
  $line=<PASSWD>;
}
close(PASSWD);
close(PASSWDOUT);
#
# get daily host times
#
$hlistfile=$logdir . "h" . $date . ".sum";
open(HLIST,"<$hlistfile");
@hlist=<HLIST>;
close(HLIST);
$dh_totaltime=0.0;
foreach $line (@hlist){
  ($host,$time)=split(' ',$line);
  $dh_time{$host}=$time;
  $dh_totaltime+=$time;
}
#
# get daily user times
#
$ulistfile=$logdir . "u" . $date . ".sum";
open(ULIST,"<$ulistfile");
@ulist=<ULIST>;
$du_totaltime=0.0;
foreach $line (@ulist){
  ($user,$time)=split(' ',$line);
  $du_time{$user}=$time;
  $du_totaltime += $time;
}
#
# get monthly host times
#
$mhlistfile=$logdir . "mh" . $month . ".sum";
open(MHLIST,"<$mhlistfile");
@mhlist=<MHLIST>;
close(MHLIST);
$header=0;
$mh_totaltime=0.0;
foreach $line (@mhlist){
  ($host,$time)=split(' ',$line);
  $mh_time{$host}=$time;
  $mh_totaltime += $time;
}
@sorted_hosts = sort { $mh_time{$b} <=> $mh_time{$a} } keys %mh_time;
#
# get monthly user times
#
$mulistfile=$logdir . "mu" . $month . ".sum";
open(MULIST,"<$mulistfile");
@mulist=<MULIST>;
close(MULIST);
$mu_totaltime=0.0;
foreach $line (@mulist){
  ($user,$time)=split(' ',$line);
  $mu_time{$user}=$time;
  $mu_totaltime += $time;
}
@sorted_users = sort { $mu_time{$b} <=> $mu_time{$a} } keys %mu_time;
#
# estimate projected usage for today and this month
#
$dh_projected_time=$dh_totaltime+($dh_totaltime/$dec_hour)*(1.0-$dec_hour);
$dh_projected_time=10*int(($dh_projected_time+5)/10);
$mh_projected_time=$mh_totaltime+($mh_totaltime/$dec_date)*(1.0-$dec_date);
$mh_projected_time=100*int(($mh_projected_time+50)/100);

$summaryfile="/var/www/html/load/summary.html";
open(SUMMARY,">$summaryfile");
print SUMMARY <<EOF;
<meta http-equiv="refresh" content="300; url=http://acrux.cfr.nist.gov/load/summary.html">
<html>
<head>
<style type="text/css">
table {
   font-family: Arial, Helvetica, sans-serif;
   font-size: 50%;
}
</style> 

<TITLE>Fire Cluster Status</title>
</head>
<BODY BGCOLOR="#FFFFFF" background="bfrlback.gif">
<small>
<strong>$date2, $ctime,
$degF</strong><br>
<a target="_blank" href="http://acrux.cfr.nist.gov/load">
Load page
</a>&nbsp;&nbsp
<a target="_blank" href="http://code.google.com/p/fds-smv/issues/list">
Issue tracker
</a>
EOF

# Temperature of $degC ($degF) measured on $clustertemptime<br>

print <<EOF;
<meta http-equiv="refresh" content="300; url=http://acrux.cfr.nist.gov/load">
<html>
<head>

<TITLE>Linux Farm Status</title>
</HEAD>
<BODY BGCOLOR="#FFFFFF" background="bfrlback.gif">

<h3>
<IMG SRC="bfrl.gif" ALT="BFRL" ALIGN="right" width="99" height="99">
Fire Research Division Linux Cluster Status:<br>
$date2, $ctime<br>
Temperature: $degF 
EOF
print <<EOF;
</h3>
<hr size="1">
<ul>
<li><a href="#currload">Current Load</a>
<li><a href="#curr">Current Usage</a>
<li><a href="#cum">Cumulative Usage</a>
<li><a href="#disk">System Status</a>
<li><a href="../../lahey">Lahey Manuals</a>
<li><a href="../../intel">Intel Manuals</a>
<li><a href="putty-0.58-installer.exe">Putty 0.58 self-installer</a>
<li><a href="faq.html">Changing passwords, adding a new user and other Cluster FAQs</a>
<li><a href="http://129.6.160.15">Temperature Sensors</a>
<li><a href="http://acrux.cfr.nist.gov/load/temp_history.png">Cluster temperaturer (last 24 hours)</a>

</ul>

<a name="currload">
<h3>Current Load</h3>
EOF
print <<EOF;
<p><table border=on>
<tr>
EOF
print SUMMARY<<EOF;
<table border=on>
<tr>
EOF
#
# determine/display unavailable systems
#
$newsystemsdown=0;
$systemcount=0;
$downsystemcount=0;
foreach $host (@hosts){
  if($systemstatus{"$host"} == 0){next;}
  $downsystemcount++;
  $newsystemsdown++;
  if($downsystemcount%7==1){
  printf <<EOF;
    </tr><tr><th>Unavailable Systems:</th>
EOF
  printf SUMMARY <<EOF;
    </tr><tr>
EOF
}
  printf <<EOF; 
<td bgcolor="black"><font size=-1 color="white">$host/</font></td>
EOF
  printf SUMMARY <<EOF; 
<td bgcolor="black"><font size=-1 color="white">$host/</font></td>
EOF
}
if($downsystemcount!=0){
printf <<EOF;
</tr>
EOF
printf SUMMARY <<EOF;
</tr>
EOF
}
open(PIDDATA,">$statdir/piddata");
open(DISKDATA,">$statdir/diskdata");
open(HOSTDATA,">$statdir/hostdata");
foreach $host (@hosts){
  if($systemstatus{"$host"}!= 0){next;}
  open(GETHOSTDATA,"</home/stats/gethdatad_$host");
  $line=<GETHOSTDATA>;
  print HOSTDATA ($line);
  $line=<GETHOSTDATA>;
  while($line ne ""){
    ($type)=split(';',$line); 
    if($type eq "pidinfo"){print PIDDATA ($line);}
    if($type eq "diskinfo"){
      ($dummy,$host,$used,$available,$usage,$disk,$disktotal)=split(';',$line);
      $hostname=$hostnames{$host};
      $newline = "$dummy;$hostname;$host;$used;$available;$usage;";
      $newline = $newline . "$disk;$disktotal";
      print DISKDATA ("$newline");
     # print DISKDATA ($line);
    }
    $line=<GETHOSTDATA>;
  }
  close (GETHOSTDATA);
}
close(PIDDATA);
close(DISKDATA);
close(HOSTDATA);

open(HOSTDATA,"<$statdir/hostdata");
$line=<HOSTDATA>;
$lastspeed="999";
$loadtotal=0.0;
$uptimesectotal=0.0;
$this_speed_total=0;
$this_speed_total2=0;
while($line ne ""){
  ($host,$upstring,$swapmem,$load1,$uptimesec,$freememory)=split(';',$line);
  $color="yellow";
  $textcolor="black";
  $uptime{$host}=$upstring;
  $freehost{$host}=$freememory;
  $swaphost{$host}=$swapmem;
  if($load1<0.5){$color="aqua";}
  if($load1>1.5){$color="red";}
#  $load1 = $load1 . " / " . 2;
  $systemload{$host}=$load1;
  $uptimesectotal+=$uptimesec;
  $loadtotal+=$load1;
  $loadcolor{$host}=$color;
  $systemcount++;
  if("$lastspeed" != "$cpuspeed{$host}"||$this_speed_total==7){
  $this_speed_total=0;
  printf <<EOF;
    </tr><tr><th>$sf $cpuspeed{$host} GHZ:$sfc</th>
EOF
}
  if($this_speed_total2==7){
  $this_speed_total2=0;
  printf SUMMARY <<EOF;
    </tr><tr>
EOF
  }
  $this_speed_total++;
  $this_speed_total2++;
  $lastspeed = $cpuspeed{$host};
  printf <<EOF; 
<td bgcolor=$loadcolor{$host}><font color=$textcolor>$sf $host - $systemload{$host}$sfc</font></td>
EOF
  printf SUMMARY <<EOF; 
<td bgcolor=$loadcolor{$host}><font color=$textcolor>$sf $host $sfc </font></td>
EOF
  $line=<HOSTDATA>;
} 
close(HOSTDATA);
$upyears=$uptimesectotal/(365.25*24*3600);
$upyears=int($upyears);
$updays=$uptimesectotal-$upyears*365.25*24*3600;
$updays=int($updays/(24*3600));
$uphours=$uptimesectotal-$upyears*365.25*24*3600-$updays*24*3600;
$uphours=int($uphours/3600);
print <<EOF;
</tr>
</table>
EOF
print SUMMARY <<EOF;
</tr>
</table>
EOF
print <<EOF;
<h3>Projected Usage</h3>
Today: $dh_projected_time h<br>
This month: $mh_projected_time h
<a name="curr">
<h3>Current Usage</h3>
EOF
#
# print process table
#
$header=0;
open(PIDDATA,"<$statdir/piddata");
@pidlines=<PIDDATA>;
close(PIDDATA);

foreach $line (sort @pidlines){
  ($dummy,$host,$user,$cputime,$exe,$pid,$memused)=split(';',$line);
  if($exe eq "lamd"){next;}
  $numhosts3{$host}++;
}
foreach $line (sort @pidlines){
  ($dummy,$host,$user,$cputime,$exe,$pid,$memused)=split(';',$line);
  if($exe eq "lamd"){next;}
  $numhosts4{$host}++;
if($header eq 0){
printf <<EOF;
<table border=on><tr>
<th>$sf System$sfc</th>
<th>$sf User$sfc</th>
<th>$sf CPU Time$sfc</th>
<th>$sf Job$sfc</th>
<th>$sf PID$sfc</th>
<th>$sf Memory$sfc</th>
</tr>
EOF
$header=1;
}
if($numhosts4{$host}==1){
printf <<EOF;
<tr><td rowspan=$numhosts3{$host}>$sf$host$sfc</td>
EOF
}
printf <<EOF;
<td>$sf$user$sfc</td>
<td align=right>$sf$cputime$sfc</td>
<td>$sf$exe$sfc</td>
<td align=right>$sf$pid$sfc</td>
<td align=right>$sf$memused$sfc</td></tr>
EOF
}
print <<EOF;
</table>
<p>
<h3>Cumulative Usage</h3>
<table>
<tr>
<td valign=top>
<table border=on>
<tr><th rowspan=2>System</th><th colspan=2>CPU Time (h)</th></tr>
<tr><th>$date3</th><th>$txtmonth</th><th></th></tr>
EOF
foreach $host (@sorted_hosts){
if($dh_time{$host}==""){
$dh_time{$host}=0.0;
}
print <<EOF;
<tr><td>$host</td><td>$dh_time{$host}</td><td>$mh_time{$host}</td></tr>
EOF
}
print <<EOF;
<tr><th>Total</th><td>$dh_totaltime</td><td>$mh_totaltime</td></tr>
<tr><th>Projected</th>
<td>$dh_projected_time</td>
<td>$mh_projected_time</td></tr>
</table>
</td>
<td valign=top>
<table border=on>
<tr><th rowspan=2>User</th><th colspan=2>CPU Time (h)</th></tr>
<tr><th>$date3</th><th>$txtmonth</th><th></th></tr>
EOF
foreach $user (@sorted_users){
if($du_time{$host}==""){
$du_time{$host}=0.0;
}
print <<EOF;
<tr><td>$user</td><td>$du_time{$user}</td><td>$mu_time{$user}</td></tr>
EOF
}
print <<EOF;
<tr><th>Total</th><td>$du_totaltime</td><td>$mu_totaltime</td></tr>
</table>
</td>
</tr>
</table>
<a name="disk">
<h3>System Status</h3>
<table border=on>
<tr>
<th rowspan=2>$sf System$sfc</th>
<th rowspan=2>$sf Load$sfc</th>
<th colspan=4>$sf Disk Usage$sfc</th>
<th colspan=2>$sf Memory (MB)$sfc</th>
<th colspan=1>$sf Swap (MB)$sfc</th>
<th rowspan=2>$sf Up Time$sfc</th>
</tr>
<tr>
<th>$sf Total<br>MB$sfc</th>
<th>$sf Free<br>MB$sfc</th>
<th>$sf Use<br>%$sfc</th>
<th>$sf Mounted On$sfc</th>
<th>$sf Total$sfc</th>
<th>$sf Free$sfc</th>
<th>$sf Free$sfc</th>
</tr>
EOF
$usedtotal=0.0;
$availabletotal=0.0;
$totaltotal=0.0;
open(DISKDATA,"<$statdir/diskdata");
@disklines=<DISKDATA>;
close(DISKDATA);
foreach $line (sort @disklines){
  ($dummy,$hostname,$host,$used,$available,$usage,$disk,$disktotal)=split(';',$line);
  $numhosts{$host}++;
}
foreach $line (sort @disklines){
  ($dummy,$hostname,$host,$used,$available,$usage,$disk,$disktotal)=split(';',$line);
  $numhosts2{$host}++;
  $availabletotal+=$available;
  $usedtotal+=$used;
  $totaltotal+=$disktotal;
  $fullfile=$logdir . "FULL" . $host;
  if($usage>90 && ! -e $fullfile ){
    system("touch $fullfile");
    open(MAIL,"|mail -s \"WARNING Disk full on $host\" $mailto");
    print MAIL <<EOF;
$date
$time
Warning the file system $disk on $host is more than 90% full.
EOF
close(MAIL);
}
if($numhosts2{$host}==1){
$hostnamexx=$hostname."(".$cpuspeed{$host}." GHZ)";
print <<EOF;
  <tr>
  <td align=right rowspan=$numhosts{$host}>$sf$hostnamexx$sfc</td>
  <td rowspan=$numhosts{$host} bgcolor=\"$loadcolor{$host}\">$sf $systemload{$host}$sfc</td>
EOF
}
print <<EOF;
  <td align=right>$sf$disktotal$sfc</td>
EOF
if($available<5000){
  $color="red";
}
else{
  if($available>20000){
    $color="aqua";
  }
else{
  $color="yellow";
}
}
if($freehost{$host}<500){
$freecolor="red";
}
else{
$freecolor="aqua";
}
printf <<EOF;
  <td align=right bgcolor=\"$color\">$sf$available$sfc</td>
  <td aligh=right bgcolor=\"$color\">$sf$usage$sfc</td>
  <td>$sf$disk$sfc</td>
EOF
if($numhosts2{$host}==1){
print <<EOF;
  <td rowspan=$numhosts{$host}>$sf $hostmem{$host}$sfc</td>
  <td bgcolor=\"$freecolor\" rowspan=$numhosts{$host}>$sf $freehost{$host}$sfc</td>
  <td rowspan=$numhosts{$host},align="right">$sf $swaphost{$host}$sfc</td>
<td rowspan=$numhosts{$host}>$sf $uptime{$host}$sfc</td>
EOF
}
printf <<EOF;
</tr>
EOF
}
$totalusage=0;
$denom=$availabletotal+$usedtotal;
if($denom!=0){
  $totalusage=$usedtotal/($availabletotal+$usedtotal);
}
$totalusage=int(100*$totalusage);
$totaltotal=int($totaltotal/1024);
$totaltotal=$totaltotal." GB";
$availabletotal=int($availabletotal/1024);
$availabletotal=$availabletotal." GB";
print <<EOF;
  <tr><td>$sf Total$sfc</td>
<td>$sf$loadtotal$sfc</td>
<td align=right>$sf$totaltotal$sfc</td>
<td align=right>$sf$availabletotal$sfc</td>
<td aligh=right>$sf$totalusage$sfc</td>
<td></td><td></td><td></td><td></td>
<td>$upyears$yearl$updays$dayl$uphours$hourl</td>
</tr>
</table>
EOF
printf <<EOF;
<p><hr size="1">
<p><a href="/">FRD Linux Farm Home</a>
</body>
</html>
EOF
printf SUMMARY <<EOF;
</body>
</html>
EOF
system("rm $lockfile");
if(-e $notifylockfile){
system("rm $notifylockfile");
}
