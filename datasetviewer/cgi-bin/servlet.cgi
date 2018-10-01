#!/usr/bin/perl
{
	use CGI;
    	use CGI::Carp qw(fatalsToBrowser);
    	use DBI;
    	use JSON;
	require './lib/UTIL.pm';
    	require './lib/AUTH.pm';
    	require './lib/Global.pm';    




    	our %VHASH = ();
    	our %CHASH = ();
    	our %SHASH = ();
    	our %GHASH = ();
    	our %SECID2NAME = ();
	our %GSECID2NAME = ();

    	$CGI::POST_MAX = 1024 * 5000;
    	our $CGI_OBJ = CGI->new();


	our $configJson;
	our $GPHASH,$PHASH;
	our $SERVER;

	$configJson->{global} = LoadConfig("conf/global.config.json");
	$configJson->{local} = LoadConfig("conf/config.json");
	$GPHASH = $configJson->{global};
	$PHASH = $configJson->{local};
	$SERVER =  $configJson->{local}{server};
       



	LoadVariables();
	$GPHASH->{device} = GetDeviceType();
        

    	print $CGI_OBJ->header();
    	if($CGI_OBJ->param('log_password')){
       		$CGI_OBJ->param(-name=>'log_password', -value=> $AUTH->_encpw($CGI_OBJ->param('log_password')));
    	}


    	$VHASH{GPAGEID} = ($VHASH{GPAGEID} ? $VHASH{GPAGEID} : $PHASH->{firstgpageid});
    	$VHASH{PAGEID} = ($VHASH{PAGEID} ? $VHASH{PAGEID} : $PHASH->{firstpageid});
        
	$VHASH{MODE} = (!$VHASH{MODE}  ? "html" : $VHASH{MODE});

    	if($VHASH{MODE} eq "json"){ #Construct and return json text
       		if($VHASH{SVC} eq "getPreviewRecords"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/getPreviewRecords.py -j '$VHASH{INJSON}'};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
                        exit;
                }
		elsif($VHASH{SVC} eq "getComment"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/getComment.py -j '$VHASH{INJSON}'};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
			exit;
		}
		elsif($VHASH{SVC} eq "saveComment"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/saveComment.py -j '$VHASH{INJSON}'};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
                        exit;
                }
		elsif($VHASH{SVC} eq "getDataModelTable"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/getDataModelTable.py};
			my $jsonText = `$cmd`;
			print "$jsonText";
			exit;
		}
		elsif($VHASH{SVC} eq "get_objects"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_objects.py};
			my $jsonText = `$cmd`;
			#print $cmd;
			print "$jsonText";
			exit;
		}
		elsif($VHASH{SVC} eq "get_single_object"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_single_object.py -o $VHASH{OBJID}};
			my $jsonText = `$cmd`;
			print "$jsonText";
			exit;
		}	
		elsif($VHASH{SVC} eq "save_object"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/save_object.py -j '$VHASH{INJSON}'};
			my $jsonText = `$cmd`;
			print "$jsonText";
			exit;
		}
                elsif($VHASH{SVC} eq "get_readme_txt"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_readme_txt.py -o $VHASH{OBJID}};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
                        exit;
                }
                elsif($VHASH{SVC} eq "get_dataset"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_dataset.py -o $VHASH{OBJID}};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
                        exit;
                }
                elsif($VHASH{SVC} eq "get_module_version"){
                    my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_module_version.py};
                    my $jsonText = `$cmd`;
                    print "$jsonText";
                    exit;
                }
        }




	my %jsParHash = (
                "deviceType"=> "device"
                ,"moduleName" => "module"
                ,"moduleBase" => "baseurl"
                ,"ghtmlRoot" => "ghtmlroot"
                ,"htmlRoot" => "htmlroot"
                ,"jsonPath" => "jsonpath"
		,"svcPath" => "svcpath"
                ,"cgiRoot" => "cgiroot"
                ,"baseUrl" => "baseurl"
        	,"moduleMenuBg" => "modulemenubg"
                ,"moduleMenuFg" => "modulemenufg"
		,"pageBase" => "pagebase"
                ,"downloadBase" => "downloadbase"
                ,"moduleRelease" => "release"
		,"server" => "server"
	);
        my %jsVarHash = (
                "gpageId" => "GPAGEID"
                ,"pageId" => "PAGEID"
		,"objId" => "OBJID"
		,"readOnly" => "READONLY"
        );
	
	my @gjsFiles = ('jquery.min.js', 'loader.js', 'vjGoogleChart.js', 'common.js', 
			 'libFuncs.js'); 
        my @gcssFiles = ('global.css', 'googlefonts.css');
        my @jsFiles = ('module.js');
        my @cssFiles = ();

	my $gheadLinks = GetGlobalHeadLinks(\@gcssFiles, \@gjsFiles);
	my $headLinks = GetModuleHeadLinks(\@cssFiles, \@jsFiles, \%jsParHash, \%jsVarHash);
	#my $headerDivOne = getHeaderDivOne();        
        my $headerDivTwo = getHeaderDivTwoNew(); 
	my $glygenLinks = getGlygenLinks();
	
	my $glygenDomain = qq{http://$SERVER.glygen.org};
	$glygenDomain = ($SERVER eq "prd" ? qq{http://glygen.org} : $glygenDomain);


	print qq{<html>
    		<head>
		<meta charset="utf-8">
    		<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
	        
    		<meta http-equiv="X-UA-Compatible" content="IE=edge">
    		<META HTTP-EQUIV="Pragma" CONTENT="no-cache">

    		<link rel="stylesheet" type="text/css" href="$glygenDomain/libraries/bootstrap/css/bootstrap.min.css">
    		<!--[if lt IE 9]>
    			<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
    			<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
    		<![endif]-->
    		<link rel="stylesheet" type="text/css" href="$glygenDomain/css/style.css">


		$gheadLinks $headLinks $glygenLinks
		
		</head> 
	 	<BODY style="background:#fff;">
		
		<script src="$glygenDomain/libraries/w3.js"></script>
    		<div w3-include-html="/content/header.html"></div>
    		<script>w3.includeHTML();</script>

		<DIV class="container" style="margin-bottom: 480px;">
			<form name=form1 method=POST action=biomuta enctype="multipart/form-data">
			$headerDivTwo
			<div class=pagewrapper>
                		<div class=pagecn id=pagecn>
				</div>
        		</div>
    			<input type=hidden name=gpageid value="$VHASH{GPAGEID}">
        		<input type=hidden name=action value="">
        		</form>
 		</DIV>
		
    		<div w3-include-html="/content/footer.html"></div>
    		<script type="text/javascript" src="$glygenDomain/libraries/jquery/jquery-3.1.js"></script>
    		<script type="text/javascript" src="$glygenDomain/libraries/bootstrap/js/bootstrap.min.js"></script>
    		<script src="$glygenDomain/js/navbar.js"></script>
    		<script src="$glygenDomain/js/ws_url.js"></script>    
		<!-- script for getting Web Service URLs -->
    		<script src="$glygenDomain/js/activity_tracker.js"></script>  
		<!-- script for logging activity -->
    		<script src="$glygenDomain/js/utility.js"></script>

		</BODY>
    		</html>
    	};


    	exit;
}














#############################
sub getGlygenLinks{


return qq{

};



}



###########################
sub getGlygenHeader{

	my $s1 =qq{width:100%;height:115px;border-top:5px solid #555;box-shadow: 2px 2px 2px #ccc;};
	my $s2 =qq{position:absolute;left:15;top:15;font-size:35;color:#333;font-family:'Arial Narrow', Arial, sans-serif;};

	my $s3 = qq{position:absolute;left:15;top:75;width:40;height:3;background:#333;};
	my $s4 = qq{position:absolute;left:15;top:90;font-size:12px;font-style:italic;color:#777;};

return qq{
<div style="$s1">
<a class="" href="http://tst.glygen.org">
<div style="$s2">GlyGen</div>
</a>
<div style="$s3"></div>
<div style="$s4">Integration of biomedical databases</div>

</div>
  




};



}
