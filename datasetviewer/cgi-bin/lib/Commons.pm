#!/usr/bin/perl


###########################
sub LoadGlobalContents{

    my @lines = `ls $PHASH{HTML_PATH}/content/page.glbl.*.html`;
    foreach my $line (@lines){
       chomp($line);
       my @parts = split(/\//, $line);
       my $fileName = $parts[-1];
       $CHASH{$fileName} = `cat $PHASH{HTML_PATH}/content/$fileName`;
       foreach my $param (keys%PHASH){
          $CHASH{$fileName} =~ s/$param/$PHASH{$param}/g;
       }
   }

}



###############################
sub GetHeaderRowOne{


    my $googleTracking = qq{<script>
        (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
        (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
        m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
        })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
        ga('create', 'UA-98197347-1', 'auto');
        ga('send', 'pageview');
        </script>};


    my $style1 = qq{position:relative;width:500px;height:60px;border:0px solid;};
    my $style2 = qq{position:relative;width:500px;height:60px;border:0px solid;};    
    my $s1 = qq{position:absolute;width:100px;height:60px;left:25px;top:20px;border:0px solid;};
    my $s2 = qq{position:absolute;width:100px;height:40px;left:110px;top:20px;border:0px solid;};
    my $cn = qq{<DIV style="$style1">
              <DIV style="$style2">
                <DIV style="$s1"> $googleTracking
		<img src="$PHASH{GHTML_ROOT}/imglib/gwu.png" height=100%></DIV>
                <DIV style="$s2"><img src="$PHASH{GHTML_ROOT}/imglib/smhs.png"></DIV>
                </DIV>
              </DIV>};
   return qq{<tr height=100 class=header><td class=header colspan=2 valign=top>$cn</td></tr>};
}


################################
sub GetHeaderRowTwo{
    

   my $hiveTitle = qq{<u>H</u>igh-performance <u>I</u>ntegrated <u>V</u>irtual <u>E</u>nvironment};
   my $cn = qq{<a href="https://hive.biochemistry.gwu.edu/dna.cgi?cmd=main" style="text-decoration:none;">
	<table cellpadding=0 style="background:#f1f1f1;" cellspacing=0 border=0>
                      <tr><td width=60>
                       <img src="$PHASH{GHTML_ROOT}/imglib/hivelogo.2.png" width=50>
                       </td><td style="color:#d0b58b;font-weight:bold;font-size:13px;">$hiveTitle</td></tr>
                      </table>
		</a>};
  return qq{<tr height=35><td colspan=2 style="background:#f1f1f1;padding:0 0 0 40;border-bottom:1px solid #ccc;" >$cn</td></tr>};

}


###############################
sub GetHeaderRowOneMobile{


    my $style1 = qq{position:relative;width:500px;height:45px;border:0px solid;};
    my $style2 = qq{position:relative;width:500px;height:45px;border:0px solid;};
    my $s1 = qq{position:absolute;width:100px;height:30px;left:5px;top:8px;border:0px solid;};
    my $s2 = qq{position:absolute;width:100px;height:30px;left:50px;top:8px;border:0px solid;};
    my $cn = qq{<DIV style="$style1">
              <DIV style="$style2">
                <DIV style="$s1"><img src="$PHASH{GHTML_ROOT}/imglib/gwu.png" height=100%></DIV>
                <DIV style="$s2"><img src="$PHASH{GHTML_ROOT}/imglib/smhs.png" height=100%></DIV>
                </DIV>
              </DIV>};
   my $row1 = qq{<tr height=20><td class=header colspan=2 valign=top>$cn</td></tr>};

   my $hiveTitle = qq{<u>H</u>igh-performance <u>I</u>ntegrated <u>V</u>irtual <u>E</u>nvironment};
   my $cn = qq{<a href="https://hive.biochemistry.gwu.edu/dna.cgi?cmd=main" style="text-decoration:none;">
        <table cellpadding=0 style="background:#f1f1f1;" cellspacing=0 border=0>
                      <tr><td width=40>
                       <img src="$PHASH{GHTML_ROOT}/imglib/hivelogo.2.png" width=35>
                       </td><td style="color:#d0b58b;font-weight:bold;font-size:11px;">$hiveTitle</td></tr>
                      </table>
                </a>};
  my $row2 = qq{<tr height=25><td colspan=2 style="background:#f1f1f1;padding:0 0 0 5;border-bottom:1px solid #ccc;" >$cn</td></tr>};

  return qq{$row1 $row2};

}




############################
sub GetSectionsGrouped{

    my @parts = split(/\;/, $PHASH{SECTIONS});
    my $rows = "";
    foreach my $part (@parts){
       my ($title, @sections) = split(/[\:,\,]/, $part);
       my $innerrows = qq{<tr><td class=sectiontitle colspan=2>$title</td></tr>};
       foreach my $comboid (@sections){
          my ($secName, $secId, $action,$access) = split(/\^/, $comboid);
          my $onclick = "form1.pageid.value = '$secId';";
          $onclick .= "form1.action.value ='$action';";
          my $class = ($secId eq $VHASH{PAGEID} ? "section2" : "section1");
          my $args = qq{pageid=$secId&action=$action};
          my $link = qq{<a href="?$args" class=$class>$secName</a>};
          my $row = qq{<tr><td width=15>&nbsp;</td><td>$link</td></tr>};
          $innerrows .= $row;
       }
       $innerrows .= qq{<tr height=20><td colspan=2>&nbsp;</td></tr>};
       $rows .=  $innerrows;
    }

    return qq{<table border=0>$rows</table>};
}   


################################
sub GetSections{
    
    my @parts = split(/\;/, $PHASH{SECTIONS});
    my $rows = "";
    foreach my $part (@parts){
       my ($title, @sections) = split(/[\:,\,]/, $part);
       my $innerrows = qq{};
       foreach my $comboid (@sections){
          my ($secName, $secId, $action,$access) = split(/\^/, $comboid);
          my $onclick = "form1.pageid.value = '$secId';";
          $onclick .= "form1.action.value ='$action';";
          my $class = ($secId eq $VHASH{PAGEID} ? "section2" : "section1");
          my $args = qq{gpageid=$VHASH{GPAGEID}&pageid=$secId};
          $args .= ($action ? "&action=$action" : "");
	  my $link = qq{<a href="?$args" class=$class id="page$secId">$secName</a>};
          my $row = qq{<tr height=20><td>$link</td></tr>};
          $innerrows .= $row;
	  $SECID2NAME{$secId} = $secName;
	  $SECID2ACTION{$secId} = $action;
          $SECID2ACCESS{$secId} = $access;
	}
       $rows .=  $innerrows;
    }
 
    return qq{<table border=0>$rows</table>};
}


################################
sub GetSectionsJson{

    my @parts = split(/\;/, $PHASH{SECTIONS});
    my $json = "[";
    my @list = ();
    foreach my $part (@parts){
       my ($title, @sections) = split(/[\:,\,]/, $part);
       my $innerrows = qq{};
       foreach my $comboid (@sections){
          my ($secName, $secId, $action,$access) = split(/\^/, $comboid);
          push @list, qq{ {"secId":"$secId",  "secName":"$secName", "action":"$action", "access":"$access",}};
          $SECID2NAME{$secId} = $secName;
          $SECID2ACTION{$secId} = $action;
       	  $SECID2ACCESS{$secId} = $access;
	}
    }
    my $json = "[" . join(",", @list) . "]";
    return $json;
}



####################
sub GetGlobalSectionsMobile{

    my @parts = split(/\;/, $PHASH{GSECTIONS});
    my $rows = qq{};
    foreach my $part (@parts){
       my ($title, @sections) = split(/[\:,\,]/, $part);
       foreach my $comboid (@sections){
          my ($secName, $secId, $action,$access) = split(/\^/, $comboid);
          my $onclick = "form1.pageid.value = '$secId';";
          $onclick .= "form1.action.value ='$action';";
          my $class = ($secId eq $VHASH{GPAGEID} ? "gsectionmobile2" : "gsectionmobile1");
          #my $url = qq{$PHASH{GHOME_URL}?gpageid=$secId&action=$action};
          my $url = ($PHASH{VER} eq "tst" ? qq{tst/$secId} : qq{$secId});
          my $link = qq{<a href="/$url" class=$class>$secName</a>};
          $rows .= ($access ne "000" ? qq{<tr height=30><td style="padding:0 0 0 20px;">$link</td></tr>} : "");
          $GSECID2NAME{$secId} = $secName;
          $GSECID2ACTION{$secId} = $action;
       }
   }
   my $secTable = qq{<table cellpadding=0 cellspacing=0 border=0>
			<tr height=50><td></td></tr>$rows</table>};

   my $prefix = ($PHASH{VER} eq "tst" ? "/tst" : "");
   my $loginlink = qq{<a href="$prefix/login" class=regular13>Login</a>};
   my $logoutlink = qq{<a href="$prefix/logout" class=regular13>Logout</a>};
   my $helplink = qq{<a href="$prefix/sitehelp" class=regular13>Help</a>};
   my $contactuslink = qq{<a href="$prefix/contactus" class=regular13>Contact Us</a>};
   my $registerlink = qq{<a href="$prefix/register" class=regular13>Register</a>};
   
   #Logout if fullname is null for some odd reason
   if(!$VHASH{FULLNAME}){
     $AUTH->logout();
   }
   my $loginStatus = $AUTH->loggedIn();
   my $msg = ($loginStatus ? qq{$VHASH{FULLNAME} &nbsp;|&nbsp; $logoutlink &nbsp;|&nbsp; $helplink &nbsp;|&nbsp; $contactuslink} :
                qq{$loginlink &nbsp;|&nbsp; $registerlink &nbsp;|&nbsp; $helplink &nbsp;|&nbsp; $contactuslink});

    #my $style = qq{position:fixed;left:0px;top:95px;width:200px;height:500px;z-index:30;display:none;};
    my $style = qq{position:absolute;left:0px;top:95px;width:200px;height:500px;z-index:30;display:none;};
    $style .= qq{background:#005581;opacity:0.95;};
    return qq{<DIV id=gsectionscn style="$style">$secTable</DIV>};

}





####################
sub GetGlobalSections{


    my @parts = split(/\;/, $PHASH{GSECTIONS});
    my $cols = qq{};
    foreach my $part (@parts){
       my ($title, @sections) = split(/[\:,\,]/, $part);
       foreach my $comboid (@sections){
          my ($secName, $secId, $action,$access) = split(/\^/, $comboid);
          my $onclick = "form1.pageid.value = '$secId';";
          $onclick .= "form1.action.value ='$action';";
          my $class = ($secId eq $VHASH{GPAGEID} ? "gsection2" : "gsection1");
          #my $url = qq{$PHASH{GHOME_URL}?gpageid=$secId&action=$action};
	  my $url = ($PHASH{VER} eq "tst" ? qq{tst/$secId} : qq{$secId});
	  my $link = qq{<a href="/$url" class=$class>$secName</a>};
          $cols .= ($access ne "000" ? qq{<td style="padding:0 50 0 0;">$link</td>} : "");
          $GSECID2NAME{$secId} = $secName;
          $GSECID2ACTION{$secId} = $action;
       }
   }
   
   my $secTable = qq{<table cellpadding=0 cellspacing=0 border=0>
                      <tr>$cols</tr></table>};

   my $prefix = ($PHASH{VER} eq "tst" ? "/tst" : "");
   my $loginlink = qq{<a href="$prefix/login" class=regular13>Login</a>};
   my $logoutlink = qq{<a href="$prefix/logout" class=regular13>Logout</a>};
   my $helplink = qq{<a href="$prefix/sitehelp" class=regular13>Help</a>};
   my $contactuslink = qq{<a href="$prefix/contactus" class=regular13>Contact Us</a>};
  
   my $registerlink = qq{<a href="$prefix/register" class=regular13>Register</a>};




   #Logout if fullname is null for some odd reason
   if(!$VHASH{FULLNAME}){
     $AUTH->logout();
   }
   my $loginStatus = $AUTH->loggedIn();
   my $msg = ($loginStatus ? qq{$VHASH{FULLNAME} &nbsp;|&nbsp; $logoutlink &nbsp;|&nbsp; $helplink &nbsp;|&nbsp; $contactuslink} : 
		qq{$loginlink &nbsp;|&nbsp; $registerlink &nbsp;|&nbsp; $helplink &nbsp;|&nbsp; $contactuslink});

   my $cn = qq{<table width=100% height=100% class=nav cellspacing=0 border=0 cellpadding=0>
               <tr>
              <td style="padding:0 0 5 40;" valign=bottom>$secTable</td>
                <td align=right  style="padding:0 15 5 0;font-size:13px;" valign=bottom>$msg</td>
                </tr>
                </table>};
   return qq{<tr height=35 style="background:#f1f1f1;"><td colspan=2 valign=bottom>$cn</td></tr>};
}

########################
sub GetNavRowMobile{

   my $homeIcon = qq{};
   my $url = qq{?pageid=10};
   my $link0 = qq{<a href="$url" class=regular13 style="font-size:11px;">$homeIcon</a>};

   my $class = ($PHASH{MODULE} eq  "LabPages"  ? "navlink2" : "navlink1");
   my $url = ($PHASH{VER} eq "tst" ? qq{tst/$VHASH{GPAGEID}} : qq{$VHASH{GPAGEID}});
   my $link1 = qq{<a href="/$url" class=$class>$GSECID2NAME{$VHASH{GPAGEID}}</a>};
   $link1 = ($VHASH{GPAGEID} ? ">> $link1" : "");

   my $class = ( $VHASH{PAGEID} eq  $PHASH{FIRSTPAGEID}  ? "navlink2" : "navlink1");
   my $url = qq{/$PHASH{BASEURL}};
   my $link2 = qq{<a href="$url" class=$class>$PHASH{MODULE}</a>};
   $link2 = ($VHASH{PAGEID} ? ">> $link2" : "");

   my $class = "navlink2";
   my $url = qq{$PHASH{BASEURL}/$VHASH{PAGEID}};
   $url = ($SECID2ACCESS{$VHASH{PAGEID}} == "000" ? "#" : $url);
   my $link3 = qq{<a href="$url" class=$class>$SECID2NAME{$VHASH{PAGEID}}</a>};
   #$link3 = (($VHASH{PAGEID} and $VHASH{PAGEID} ne  $PHASH{FIRSTPAGEID}) ? ">> $link3" : "");
   $link3 = ($VHASH{PAGEID} ? ">> $link3" : "");

   my $style = qq{font-size:10px; padding:0 0 5 10;background:#f1f1f1;border-bottom:1px solid #ccc;};
   my $cn = qq{$link1 $link2 $link3};

   return qq{<tr height=25><td style="$style" valign=bottom colspan=2>$cn</td></tr>};
}



########################
sub GetNavRow{

   my $homeIcon = qq{};
   my $url = qq{?pageid=10};
   my $link0 = qq{<a href="$url" class=regular13 style="font-size:11px;">$homeIcon</a>};

   my $class = ($PHASH{MODULE} eq  "LabPages"  ? "navlink2" : "navlink1");
   my $url = ($PHASH{VER} eq "tst" ? qq{tst/$VHASH{GPAGEID}} : qq{$VHASH{GPAGEID}});
   my $link1 = qq{<a href="/$url" class=$class>$GSECID2NAME{$VHASH{GPAGEID}}</a>};
   $link1 = ($VHASH{GPAGEID} ? ">> $link1" : ""); 

   my $class = ( $VHASH{PAGEID} eq  $PHASH{FIRSTPAGEID}  ? "navlink2" : "navlink1");
   my $url = qq{$PHASH{BASEURL}};
   my $link2 = qq{<a href="$url" class=$class>$PHASH{MODULE}</a>};
   $link2 = ($VHASH{PAGEID} ? ">> $link2" : "");    


   my $class = "navlink2";
   my $url = qq{$PHASH{BASEURL}/$VHASH{PAGEID}};
   $url = ($SECID2ACCESS{$VHASH{PAGEID}} == "000" ? "#" : $url);
   my $link3 = qq{<a href="$url" class=$class>$SECID2NAME{$VHASH{PAGEID}}</a>};
   $link3 = ($VHASH{PAGEID} ? ">> $link3" : "");
  
 
   my $style = qq{font-size:10px; padding:0 0 5 40;background:#f1f1f1;border-bottom:1px solid #ccc;};
   my $cn = qq{$link1 $link2 $link3};

   return qq{<tr height=30><td style="$style" valign=bottom colspan=2>$cn</td></tr>};
}



############################
sub GetWelcomeRow{

   my $loginlink = qq{<a href="login" class=regular13>Login</a>};
   my $logoutlink = qq{<a href="logout" class=regular13>Logout</a>};
   my $msg = ($AUTH->loggedIn ? qq{Welcome $VHASH{USERNAME}! $logoutlink} : $loginlink);
  
   my $link1 = qq{<a href="$PHASH{DASHBOARDURL}" class=regular13>HIVE Lab</a>}; 
   my $cn = qq{<table width=100% height=100% class=nav cellspacing=0 border=0>
               <tr>
              <td class=nav >$link1 >> $PHASH{MODULE} >> $GSECID2NAME{$VHASH{PAGEID}}   </td>
                <td align=right class=welcome valign=bottom>$msg</td>
                </tr>
                </table>};
   return qq{<tr class=nav><td colspan=2>$cn</td></tr>};

}



##########################
sub LoadUserGroups{
 
    my $sql = qq{SELECT * FROM auth_group WHERE userid = '$VHASH{USERID}'};
    my $refs = $DBH->selectall_arrayref($sql, {Slice => {}});
    foreach my $ref (@$refs){
       $GHASH{$ref->{groupname}} = 1;
    }
}




######################################
sub SaveUploadedFile{
    my($name, $savedir, $savefilename) = @_;

    my $srcfilename = $CGI_OBJ->param($name);
    if (!$srcfilename )  {
        $VHASH{ERR} = "Please select file!<br>";
        return;
    }

    my @parts = split(/\./, $srcfilename);
    my $ext = lc(pop(@parts));
    my $newfilename = "$savedir/$savefilename.$ext";
    my $tmpfilename = $CGI_OBJ->tmpFileName($srcfilename);

    `mv  $tmpfilename $newfilename`;
    `chmod 755 $newfilename`;
    return $newfilename;
}


###########################
sub GetGlobalMenuIcon{

    my $style1 = qq{left:0px;top:0px;width:100%;height:100%;};
    my $style2 = qq{position:absolute;left:0px;width:100%;height:4px;background:#fff;};
    #$style2 .= qq{-moz-border-radius: 15px;border-radius: 15px;};

    my $cn = qq{<DIV style="$style1">};
    $cn .= qq{<DIV style="$style2;top:0px;"></DIV>};
    $cn .= qq{<DIV style="$style2;top:9px;"></DIV>};
    $cn .= qq{<DIV style="$style2;top:18px;"></DIV>};
    $cn .= qq{</DIV>};

    #my $style = qq{position:fixed;right:10px;top:12px;width:6px;height:25px;z-index:20;};
    my $style = qq{position:absolute;right:10px;top:12px;width:6px;height:25px;z-index:20;};
    return qq{<a href="" id=gmenuicon><DIV style="$style">$cn</DIV></a>};                

}


###########################
sub GetModuleIcon{
    my($radius, $hh, $yjump, $bgcolor) = @_;
 
    if(!$radius){
        $radius = 15 . "px";
	$hh = 4 . "px";
	$yjump = 7;
    }
	
    my $style1 = qq{left:0px;top:0px;width:100%;height:100%;};
    my $style2 = qq{position:absolute;left:0px;width:100%;height:$hh;background:$bgcolor;};
    $style2 .= qq{-moz-border-radius: $radius;border-radius: $radius;};

    my $cn = qq{<DIV style="$style1">};
    my $top = "0px";
    $cn .= qq{<DIV style="$style2;top:$top;"></DIV>};
    my $top = $yjump . "px";
    $cn .= qq{<DIV style="$style2;top:$top;"></DIV>};
    my $top = 2*$yjump . "px";
    $cn .= qq{<DIV style="$style2;top:$top;"></DIV>};
    $cn .= qq{</DIV>};
    
   return $cn;


}


#########################
sub LoadPageContents{
    my($htmlPath, $hashptr) = @_;

    
    foreach my $fileName (@{$hashptr}){
        $CHASH{$fileName} = `cat $htmlPath/$fileName`;
        foreach my $param (keys%PHASH){
          $CHASH{$fileName} =~ s/$param/$PHASH{$param}/g;
       }
   }
}



#############################################
sub GetHeaderDivOne(){


        my $gsections = ($PHASH{DEVICE} eq "pc" ? GetGlobalSections() : "");
        my $gsectionsMenu = ($PHASH{DEVICE} eq "pc" ? "" : GetGlobalSectionsMobile());
        my $sections = GetSections();
        my $navRow = ($PHASH{DEVICE} eq "pc" ? GetNavRow() : GetNavRowMobile());
        my $headerRowOne = "";
        my $headerRowTwo = "";
        my $gmenuIconDiv = "";
        if($PHASH{DEVICE} ne "pc"){
                $gmenuIconDiv = GetGlobalMenuIcon();
                $headerRowOne = GetHeaderRowOneMobile();
        }
        else{
                $headerRowOne = GetHeaderRowOne();
                $headerRowTwo = GetHeaderRowTwo();
        }
        my $style5 = qq{position:absolute;left:0px;width:100%;height:25px;font-size:10px;};
        $style5 .= qq{font-weight:bold;color:#777;text-align:right;z-index:9;padding:0px 0px 0px 0px;};
        $style5 .= ($PHASH{VER} eq "tst" ? "background:#F9E79F;" : "background:#FFFFFF;");
        $style5 .= ($PHASH{DEVICE} eq "pc" ? "top:200px;" : "top:95px;");
        my $versionDiv = qq{<DIV style="$style5">Version - $PHASH{VER} &nbsp;&nbsp;&nbsp;&nbsp; &nbsp; </DIV>};

        my $table1 = qq{<table width=100% cellpadding=0 cellspacing=0 border=0>
                        $headerRowOne $headerRowTwo $gsections $navRow
                        </table>};
        my $div = qq{<DIV class=header style="">$table1 $versionDiv</DIV>};

        return $div;
}


#############################################
sub GetHeaderDivTwo(){
	my($moduleMenuBg, $moduleMenuFg) = @_;


        my $style3 = qq{position:absolute;width:20px;height:25px;z-index:11;};
        $style3 .= ($PHASH{DEVICE} eq "pc" ? "top:23px;" : "top:75px;");
        $style3 .= ($PHASH{DEVICE} eq "pc" ? "left:15px;" : "left:10px;");

        my $style4 = qq{position:absolute;left:0px;width:250px;height:500px;};
        $style4 .= qq{padding:0px 0px 0px 0px;background:$moduleMenuBg;display:none;border:0px solid red;z-index:10;};
        $style4 .= ($PHASH{DEVICE} eq "pc" ? "top:50px;" : "top:95px;");

        my $icon = GetModuleIcon(15, 3, 5, $moduleMenuFg);
        my $moduleIcondiv = qq{<a href="" id=moduleicon><DIV style="$style3">$icon</DIV></a>};
        my $moduleSectionsdiv = qq{<DIV style="$style4" id=sectionscn></DIV>};
        my $cellstyl1 = qq{padding:20 0 10 50;background:$moduleMenuBg;color:$moduleMenuFg;};
        my $tableA = qq{<table  width=100% cellpadding=0 cellspacing=0 border=0>
                        <tr height=35><td class=logotext1 style="$cellstyl1">$PHASH{TITLE}</td></tr>
                        </table>};
	my $cellstyl2 = qq{padding:10 0 0 20;background:$moduleMenuBg;color:$moduleMenuFg;font-size:12px;};
	my $tableB = qq{<table class=global width=100% cellpadding=0 cellspacing=0 border=0>
		<tr height=40><td style="$cellstyl2" id="searchboxcn" valign=top>&nbsp;</td></tr>
		<tr class=menu><td class=pagecn id="pagecn" style="padding:40px 0px 300px 0px;" colspan=3></td></tr>
		</table>};
	my $style = qq{position:absolute;height:40px;z-index:102;border:0px solid red;};
        $style .= ($PHASH{DEVICE} eq "pc" ? "top:220px;" : "top:100px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "left:50px;" : "left:0px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "width:80%;" : "width:100%;");
        my $divA = qq{<DIV style="$style">$tableA $moduleIcondiv $moduleSectionsdiv</DIV>};

	my $style = qq{position:relative;height:100%;border:0px solid red;};
        $style .= ($PHASH{DEVICE} eq "pc" ? "top:70px;" : "top:300px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "left:50px;" : "left:0px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "width:80%;" : "width:100%;");
	my $divB = qq{<DIV style="$style">$tableB</DIV>};
        return qq{$divA $divB};
}

#############################################
sub GetHeaderDivTwoNoSearchBox(){
        my($moduleMenuBg, $moduleMenuFg, $rc_html) = @_;


        my $style3 = qq{position:absolute;width:20px;height:25px;z-index:11;};
        $style3 .= ($PHASH{DEVICE} eq "pc" ? "top:23px;" : "top:75px;");
        $style3 .= ($PHASH{DEVICE} eq "pc" ? "left:15px;" : "left:10px;");

        my $style4 = qq{position:absolute;left:0px;width:250px;height:500px;};
        $style4 .= qq{padding:0px 0px 0px 0px;background:$moduleMenuBg;display:none;border:0px solid red;z-index:10;};
        $style4 .= ($PHASH{DEVICE} eq "pc" ? "top:50px;" : "top:95px;");

        my $icon = GetModuleIcon(15, 3, 5, $moduleMenuFg);
        my $moduleIcondiv = qq{<a href="" id=moduleicon><DIV style="$style3">$icon</DIV></a>};
        my $moduleSectionsdiv = qq{<DIV style="$style4" id=sectionscn></DIV>};
        my $cellstyl1 = qq{padding:0 0 10 50;background:$moduleMenuBg;color:$moduleMenuFg;};
        my $tableA = qq{<table  width=100% cellpadding=0 cellspacing=0 border=0>
                        <tr height=70><td class=logotext1 style="$cellstyl1">$PHASH{TITLE}</td></tr>
                        </table>};
        my $cellstyl2 = qq{padding:0 0 0 20;background:$moduleMenuBg;color:$moduleMenuFg;font-size:12px;};
        

	my $rcRow = ($rc_html ne "" ? qq{<tr><td style="padding:40 50 0 30;" >$rc_html</td></tr>} :"");
	my $style = ($rc_html eq "" ? "padding:40px 0px 300px 0px;" : "padding:0px 0px 300px 0px;");
	my $tableB = qq{<table class=global width=100% cellpadding=0 cellspacing=0 border=0>
                $rcRow
		<tr class=menu><td class=pagecn id="pagecn" style="$style" colspan=3></td></tr>
                </table>};

	my $style = qq{position:absolute;height:60px;z-index:102;border:0px solid red;};
        $style .= ($PHASH{DEVICE} eq "pc" ? "top:220px;" : "top:100px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "left:30px;" : "left:0px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "width:90%;" : "width:100%;");
        my $divA = qq{<DIV style="$style">$tableA $moduleIcondiv $moduleSectionsdiv</DIV>};

        my $style = qq{position:relative;height:100%;border:0px solid red;};
        $style .= ($PHASH{DEVICE} eq "pc" ? "top:90px;" : "top:300px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "left:30px;" : "left:0px;");
        $style .= ($PHASH{DEVICE} eq "pc" ? "width:90%;" : "width:100%;");
        my $divB = qq{<DIV style="$style">$tableB</DIV>};
        return qq{$divA $divB};
}





###############################################
sub GetHeadLinks(){
        my($cssFileListPtr, $jsFileListPtr, $jsParHashPtr, $jsVarHashPtr) = @_;

        my $viewPort = qq{<meta name="viewport" content ="width=device-width, initial-scale=1">};
        $viewPort = ($PHASH{DEVICE} eq "pc" ? "" : $viewPort);
        my $sectionsJson = GetSectionsJson();

        my $cssLinks = qq{};
        foreach my $fname (@$cssFileListPtr){
                $cssLinks .= qq{<link rel="stylesheet" href="$PHASH{GHTML_ROOT}/csslib/$fname" type="text/css">\n};
        }
        my $jsLinks = qq{<script language="javascript" src="$PHASH{GHTML_ROOT}/jslib/common.js"></script>\n};
        foreach my $fname (@$jsFileListPtr){
                $jsLinks .= qq{<script language="javascript" src="$PHASH{HTML_ROOT}/jslib/$fname"></script>\n};
        }

        my $jsGlobal = qq{<script>};
        my %jsParHash = %{$jsParHashPtr};
        for my $par (keys%jsParHash){
                $jsGlobal .= qq{var $par = \'$PHASH{$jsParHash{$par}}\';\n};
        }

        my %jsVarHash = %{$jsVarHashPtr};
        for my $var (keys%jsVarHash){
                $jsGlobal .= qq{var $var = \'$VHASH{$jsVarHash{$var}}\';\n};
        }
        $jsGlobal .= qq{var sectionObjs = $sectionsJson;\n};
        $jsGlobal .= qq{</script>};
        return qq{$cssLinks $jsLinks $jsGlobal $viewPort};
}


###############################################
sub GetHeadLinksOne(){
        my($cssFileListPtr, $jsFileListPtr) = @_;

        my $viewPort = qq{<meta name="viewport" content ="width=device-width, initial-scale=1">};
        $viewPort = ($PHASH{DEVICE} eq "pc" ? "" : $viewPort);
        my $cssLinks = qq{};
        foreach my $fname (@$cssFileListPtr){
                $cssLinks .= qq{<link rel="stylesheet" href="$PHASH{GHTML_ROOT}/csslib/$fname" type="text/css">\n};
        }
        my $jsLinks = qq{};
        foreach my $fname (@$jsFileListPtr){
                $jsLinks .= qq{<script language="javascript" src="$PHASH{GHTML_ROOT}/jslib/$fname"></script>\n};
        }
        return qq{$cssLinks $jsLinks $viewPort};
}



###############################################
sub GetHeadLinksTwo(){
        my($cssFileListPtr, $jsFileListPtr, $jsParHashPtr, $jsVarHashPtr) = @_;

        my $sectionsJson = GetSectionsJson();
        my $cssLinks = qq{};
        foreach my $fname (@$cssFileListPtr){
                $cssLinks .= qq{<link rel="stylesheet" href="$PHASH{HTML_ROOT}/csslib/$fname" type="text/css">\n};
        }
        my $jsLinks = qq{};
        foreach my $fname (@$jsFileListPtr){
                $jsLinks .= qq{<script language="javascript" src="$PHASH{HTML_ROOT}/jslib/$fname"></script>\n};
        }

        my $jsGlobal = qq{<script>};
        my %jsParHash = %{$jsParHashPtr};
        for my $par (keys%jsParHash){
                $jsGlobal .= qq{var $par = \'$PHASH{$jsParHash{$par}}\';\n};
        }

        my %jsVarHash = %{$jsVarHashPtr};
        for my $var (keys%jsVarHash){
                $jsGlobal .= qq{var $var = \'$VHASH{$jsVarHash{$var}}\';\n};
        }
        $jsGlobal .= qq{var sectionObjs = $sectionsJson;\n};
        $jsGlobal .= qq{</script>};
        return qq{$cssLinks $jsLinks $jsGlobal $viewPort};
}




1;

