#!/usr/bin/perl



#################################
sub getHeaderDivOne{

	my %sHash = ();
	$sHash{header} = qq{padding:16px 0px 16px 0px; border-bottom:1px solid #fff;};
	$sHash{header} .= qq{-moz-box-shadow:0 5px 5px rgba(182, 182, 182, 0.75);};
	$sHash{header} .= qq{-webkit-box-shadow: 0 5px 5px rgba(182, 182, 182, 0.75);};
	$sHash{header} .= qq{box-shadow: 0 5px 5px rgba(182, 182, 182, 0.75);};
	$sHash{header} .= qq{position:absolute;top:0px;};
        
	return qq{
		<nav class="navbar navbar-expand-lg navbar-dark bg-dark fixed-top" style="$sHash{header}">
		<div class="container" style="width:100%;margin-left:2.1%;margin-right:2.1%;border:0px solid red;">
  		<div style="position:relative;left:13px;top:0px;width:280px;height:60px;border:0px dashed #fff;">
                        <div id="logocn1" style="position:absolute;left:0px;">
                                <img src="$GPHASH->{$SERVER}{rootinfo}{ghtmlroot}/gimglib/logo.png" width=100%>
                        </div>
                        <div id="logocn2" style="position:absolute;left:92px;">
                                <img src="$GPHASH->{$SERVER}{rootinfo}{ghtmlroot}/gimglib/smhs.png" width=100%>
                        </div>
                 </div>
        	<button class="navbar-toggler" type="button" data-toggle="collapse" 
			data-target="#navbarResponsive" aria-controls="navbarResponsive" 
			aria-expanded="false" aria-label="Toggle navigation">
          	<span class="navbar-toggler-icon"></span>
        	</button>
        	<div class="collapse navbar-collapse" id="navbarResponsive" 
			style="position:relative;top:30px;border:0px dashed #fff;">
        	</div>
      	</div>
    	</nav>
 	};
}

#################################
sub getHeaderDivTwoNew{

    #my $icon = GetModuleIcon(15, 3, 7, $PHASH->{modulemenufg});
        my $icon = "";

        my $searchBoxDiv = qq{<div class="modulesearchboxcn" id=modulesearchboxcn></div>};
        my $moduleHeaderClass = qq{moduleheaderwrapper};
        if($PHASH->{queryform}{nosearch}){
                $searchBoxDiv = qq{};
                $moduleHeaderClass = qq{moduleheaderwrappernosearch};
        }
       
        my $menuTitleDiv = qq{
                <div id=pagelinkscn style="width:100%;height:20;text-align:right;"></div>
		<div class=$moduleHeaderClass style="background:$PHASH->{modulemenubg};">
                        <div class="moduleiconcn" id=moduleiconcn>$icon</div>
                        <div id="moduletitlecn" class="moduletitlecn" style="color:$PHASH->{modulemenufg};">
                        $PHASH->{title}
                        </div>
                        <div id="moduleversioncn" class="moduleversioncn" style="color:$PHASH->{modulemenufg};">
                        </div>
                        $searchBoxDiv
                        <div class="modulesectionscn" id=modulesectionscn style="background:$PHASH->{modulemenubg};"></div>
                </div>
	};
        my $moduleTitle = ($PHASH->{module} eq "global" ? "" : $menuTitleDiv);
        return qq{
                <div class="labtitlewrapper">
                        $moduleTitle
                </div>
        };
}


#################################
sub getHeaderDivTwo{

        my $icon = GetModuleIcon(15, 3, 7, $PHASH->{modulemenufg});
	my $searchBoxDiv = qq{<div class="modulesearchboxcn" id=modulesearchboxcn></div>};
	my $moduleHeaderClass = qq{moduleheaderwrapper};
	if($PHASH->{queryform}{nosearch}){
		$searchBoxDiv = qq{};
		$moduleHeaderClass = qq{moduleheaderwrappernosearch};
	}
	my $menuTitleDiv = qq{
		<div class=$moduleHeaderClass style="background:$PHASH->{modulemenubg};">
			<div class="moduleiconcn" id=moduleiconcn>$icon</div>
			<div class="moduletitlecn" style="color:$PHASH->{modulemenufg};">$PHASH->{title}</div>
			$searchBoxDiv
			<div class="modulesectionscn" id=modulesectionscn style="background:$PHASH->{modulemenubg};"></div>
		</div>
	};
	my $moduleTitle = ($PHASH->{module} eq "global" ? "" : $menuTitleDiv);
	return qq{
                <div class="labtitlewrapper">
                        <a href="$GPHASH->{$SERVER}{rootinfo}{ghomeurl}" style="text-decoration:none;">
			<div class=labtitlecn>$GPHASH->{title}</div>
                	</a>
			$moduleTitle
		</div>
        };
}





##########################
sub LoadUserGroups{
 
    my $sql = qq{SELECT * FROM auth_group WHERE userid = '$VHASH{USERID}'};
    my $refs = $DBH->selectall_arrayref($sql, {Slice => {}});
    foreach my $ref (@$refs){
       $GHASH{$ref->{groupname}} = 1;
    }
}






###############################################
sub GetGlobalHeadLinks(){
        my($cssFileListPtr, $jsFileListPtr) = @_;

        my $viewPort = qq{<meta name="viewport" content ="width=device-width, initial-scale=1">};
        $viewPort = ($GPHASH->{device} eq "pc" ? "" : $viewPort);
        my $cssLinks = qq{};
        foreach my $fname (@$cssFileListPtr){
                $cssLinks .= qq{<link rel="stylesheet" href="$GPHASH->{$SERVER}{rootinfo}{ghtmlroot}/gcsslib/$fname" type="text/css">\n};
        }
        my $jsLinks = qq{};
        foreach my $fname (@$jsFileListPtr){
                $jsLinks .= qq{<script language="javascript" src="$GPHASH->{$SERVER}{rootinfo}{ghtmlroot}/gjslib/$fname"></script>\n};
        }
        return qq{$cssLinks $jsLinks $viewPort};
}



###############################################
sub GetModuleHeadLinks(){
        my($cssFileListPtr, $jsFileListPtr, $jsParHashPtr, $jsVarHashPtr) = @_;

	my $sections = ( $PHASH->{module} eq "global" ?  $GPHASH->{sections} : $PHASH->{sections});
	my $sectionsJson = encode_json $sections;
        my $queryFormJson = encode_json $PHASH->{queryform};
	my $rootInfoJson = encode_json $PHASH->{$SERVER}{rootinfo};
	my $grootInfoJson = encode_json $GPHASH->{$SERVER}{rootinfo};

	my $topSectionsJson = encode_json $GPHASH->{topsections};

	foreach my $obj (@{$PHASH->{sections}}){
                my ($secName, $secId) = ($obj->{label}, $obj->{id});
                $SECID2NAME{$secId} = $secName;
        }

	my $cssLinks = qq{};
        foreach my $fname (@$cssFileListPtr){
                $cssLinks .= qq{<link rel="stylesheet" href="$PHASH->{$SERVER}{rootinfo}{htmlroot}/csslib/$fname" type="text/css">\n};
        }
        my $jsLinks = qq{};
        foreach my $fname (@$jsFileListPtr){
                $jsLinks .= qq{<script language="javascript" src="$PHASH->{$SERVER}{rootinfo}{htmlroot}/jslib/$fname"></script>\n};
        }

        my $jsGlobal = qq{<script>};
        my %jsParHash = %{$jsParHashPtr};
        for my $par (keys%jsParHash){
                my $val = ($GPHASH->{$jsParHash{$par}} ? $GPHASH->{$jsParHash{$par}} :
                                        $GPHASH->{$SERVER}{rootinfo}{$jsParHash{$par}});
		$val = ($val ? $val : $PHASH->{$jsParHash{$par}});
		$val = ($val ? $val : $PHASH->{$SERVER}{rootinfo}{$jsParHash{$par}});
		$jsGlobal .= qq{var $par = \'$val\';\n};
        }

        my %jsVarHash = %{$jsVarHashPtr};
        for my $var (keys%jsVarHash){
                $jsGlobal .= qq{var $var = \'$VHASH{$jsVarHash{$var}}\';\n};
        }
        $jsGlobal .= qq{var sectionObjs = $sectionsJson;\n};
	$jsGlobal .= qq{var topSectionObjs = $topSectionsJson;\n};
	$jsGlobal .= qq{var queryFormJson = $queryFormJson;\n};
	$jsGlobal .= qq{var rootInfoJson = $rootInfoJson;\n};
	$jsGlobal .= qq{var grootInfoJson = $grootInfoJson;\n};
        $jsGlobal .= qq{var moduleVersion = \'$SERVER\';\n};
	$jsGlobal .= qq{</script>};
        return qq{$cssLinks $jsLinks $jsGlobal $viewPort};
}


###############################
sub LoadConfig{
	my ($configFile) = @_;

	open FR, $configFile;
  	my $string = "";
	while ( my $line = <FR>){
		$string .= $line
	}
  	close FR;
	my $jsonObj = decode_json $string;
	return $jsonObj;
}


###############################
sub LoadHtmlFileContent{
        my ($htmlFile) = @_;

        open FR, $htmlFile;
        my $string = "";
        while ( my $line = <FR>){
                $string .= $line
        }
        close FR;
        return $string
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



###########################
sub GetFooter{


	my @links = (
		"<a class=footer href=\"http://www.campusadvisories.gwu.edu/\">Campus Advisories</a>"
		,"<a class=footer href=\"http://www.gwu.edu/legal/copyright\">Copyright</a>"
		,"<a class=footer href=\"http://www.gwu.edu/legal/privacypolicy\" target=\"_blank\">Privacy Policy</a>"
		,"<a class=footer href=\"http://www.gwu.edu/legal/termsofuse\" target=\"_blank\">Terms of Use</a>"
		,"<a class=footer href=\"http://www.gwu.edu/contact-gw\" target=\"_blank\">Contact GW</a>"
		,"<a class=footer href=\"http://www.gwu.edu/az-index\" target=\"_blank\">A - Z Index</a>"
	);

	my $table1 = qq{<table align=right style="margin:20 0 0 0;" cellpadding=10><tr>};
	foreach my $link (@links){
		$table1 .= qq{<td>$link</td>};
	}
	$table1 .= qq{</td></tr></table>};

	my $table = qq{<table width=100%>
		<tr>
		<td style="padding:30 0 0 30;"><img src="$GPHASH->{$SERVER}{rootinfo}{ghtmlroot}/gimglib/smhs-full.png" width=50%></td>
		<td>$table1</td>
		</tr>
		</table>};
	return $table;


}



1;

