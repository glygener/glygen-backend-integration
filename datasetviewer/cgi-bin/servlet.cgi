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
		if($VHASH{SVC} eq "get_comment"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_comment.py -j '$VHASH{INJSON}'};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
			exit;
		}
		elsif($VHASH{SVC} eq "save_comment"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/save_comment.py -j '$VHASH{INJSON}'};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
                        exit;
                }
		elsif($VHASH{SVC} eq "get_data_model_table"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_data_model_table.py};
			my $jsonText = `$cmd`;
			print "$jsonText";
			exit;
		}
		elsif($VHASH{SVC} eq "search_objects"){
			my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/search_objects.py -j '$VHASH{INJSON}'};
			my $jsonText = `$cmd`;
			#print $cmd;
			print "$jsonText";
			exit;
		}
                elsif($VHASH{SVC} eq "get_readme_txt"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_readme_txt.py -o $VHASH{OBJID} -v $VHASH{OBJVER}};
                        my $jsonText = `$cmd`;
                        print "$jsonText";
                        exit;
                }
                elsif($VHASH{SVC} eq "get_bco"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_bco.py -o $VHASH{OBJID} };
                        $cmd .= ($VHASH{OBJVER} ne ""  ? " -v $VHASH{OBJVER}" : "");
                        my $jsonText = `$cmd`;
                        print "$jsonText";
                        #print qq{$cmd};
                        exit;
                }
                elsif($VHASH{SVC} eq "get_dataset"){
                        my $cmd = qq{python $PHASH->{$SERVER}{pathinfo}{cgipath}/svc/get_dataset.py -o $VHASH{OBJID} };
                        $cmd .= ($VHASH{OBJVER} ne "\"\""  ? " -v $VHASH{OBJVER}" : "");
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
	        ,"bcoPrefix" => "bcoprefix"
        );
        my %jsVarHash = (
                "gpageId" => "GPAGEID"
                ,"pageId" => "PAGEID"
		,"objId" => "OBJID"
		,"readOnly" => "READONLY"
        );
        
        my @gjsFiles = ();
        my @gcssFiles = ();
        my @jsFiles = ('jquery.min.js', 'loader.js', 'vjGoogleChart.js', 'global.js', 'module.js');
        my @cssFiles = ('global.css', 'googlefonts.css');

	my $gheadLinks = GetGlobalHeadLinks(\@gcssFiles, \@gjsFiles);
	my $headLinks = GetModuleHeadLinks(\@cssFiles, \@jsFiles, \%jsParHash, \%jsVarHash);
	#my $headerDivOne = getHeaderDivOne();        
        my $headerDivTwo = getHeaderDivTwo(); 
	

        my $moduleDomain = qq{https://$SERVER.$PHASH->{project}.org};
        $moduleDomain = ($SERVER eq "prd" ? qq{https://$PHASH->{project}.org} : $moduleDomain);
        
        my $projectHeader = getProjectHeader($moduleDomain, $PHASH->{project});
        my $projectFooter = getProjectFooter($moduleDomain, $PHASH->{project});


	print qq{<html>
    		<head>
		<meta charset="utf-8">
    		<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
	        
    		<meta http-equiv="X-UA-Compatible" content="IE=edge">
    		<META HTTP-EQUIV="Pragma" CONTENT="no-cache">
		$gheadLinks $headLinks 
		
		</head> 
	 	<BODY style="background:#fff;">

                $projectHeader

                <DIV class="container-fluid text-center homepageContainer">
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
	
                $projectFooter
	

		</BODY>
    		</html>
    	};


    	exit;
}




###########################
sub getProjectHeader{
        my($moduleDomain, $project) = @_;

        if ($project eq "oncomx"){
                return qq{
                <link href="/csslib/style.css" rel="stylesheet">    
                <link href="$moduleDomain/static/themes/Imperial/lib/bootstrap/css/bootstrap.min.css" rel="stylesheet">
                        <link href="$moduleDomain/static/themes/Imperial/css/style.css?2" rel="stylesheet">
                        <link href="$moduleDomain/static/themes/Imperial/css/onco-custom.css" rel="stylesheet">
        
                        <header id="header">
                        <div class="container">
                        <nav id="nav-menu-container" style="float:left">
                        <ul class="nav-menu">
                        <li><a href="$moduleDomain#home" class="btn-get-started pull-left">Home</a></li>
                        <li><a href="$moduleDomain#about" class="btn-get-started">Explore</a></li>
                        <li><a href="$moduleDomain#subscribe" class="btn-get-started">Downloads</a></li>
                        <li><a href=“$moduleDomain#about" class="btn-get-started">Statistics</a></li>
                        <li><a href="$moduleDomain/oncodataview/" class="btn-get-started">Original Data Sources</a></li>
                        </ul>
                        </nav>
                        <nav id="nav-menu-container" style="float:right">
                        <ul class="nav-menu">
                        <li><a href="$moduleDomain/static/docs/OncoMX_global_readme_v1.0_02Oct2018.txt">Help</a></li>
                         <li><a href="$moduleDomain#team">About Us</a></li>
                        <li><a href="$moduleDomain#contact">Contact Us</a></li>
                        </ul>
                        </nav>
                        </div>
                        </header>
                };
        }
        else{
                return qq{
                <link rel="stylesheet" type="text/css" href="$moduleDomain/libraries/bootstrap/css/bootstrap.min.css">
                <!--[if lt IE 9]>
                        <script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
                        <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
                <![endif]-->
                <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.4.1/css/all.css" integrity="sha384-5sAR7xN1Nv6T6+dT2mhtzEpVJvfS3NScPQTrOxhwjIuvcA67KV2R5Jz6kr4abQsz" crossorigin="anonymous">
                <link rel="stylesheet" type="text/css" href="$moduleDomain/css/alertify.css">
                <link rel="stylesheet" type="text/css" href="$moduleDomain/css/style.css">
                <link rel="stylesheet" type="text/css" href="$moduleDomain/css/server.css">
                <link rel="stylesheet" href="$moduleDomain/css/reset.css">

                <div>
                
                    <div class="alert navbar-fixed gg-alert fade in" id="tracking_banner">
                    
                <span>Do you want <strong>GlyGen</strong> to remember your searches for your future use? This can be changed at any time in the <strong>My GlyGen</strong> section.</span>
                
                <br>
                
                <button type="button" class="btn btn-default gg-btn-margin" onclick='logID()'>Allow</button>
                   
                <button type="button" class="btn btn-default gg-btn-margin" onclick='doNotLog()'>Don't Allow</button>
                
                </div>
                <!----> 
                <nav class="navbar navbar-inverse gg-navbar gg-navbar-inverse">
                    <div class="container-fluid">
                        <div class="navbar-header"> 
                            <button type="button" class="navbar-toggle" data-toggle="collapse" data-target="$moduleDomain/index.html#myNavbar">
                                <span class="icon-bar"></span>
                                <span class="icon-bar"></span>
                                <span class="icon-bar"></span>                   
                            </button>
                            <a class="navbar-brand glygenmenu" id="index.html" href="$moduleDomain/index.html"><img id="logo" src="/content/logo-glygen.svg">
                            </a>
                        </div>
                        
                        <div class="collapse navbar-collapse gg-collapse gg-nav-item" id="myNavbar"> 
                            <ul class="nav navbar-nav gg-nav-nav ">
                                <li>
                                    <a alt="home" class=glygenmenu href="$moduleDomain/index.html">HOME</a>
                                </li>
                                <li class="dropdown">
                                    <a classs="dropdown-toggle glygenmenu" data-toggle="dropdown" alt="explore" href="$moduleDomain/index.html">EXPLORE 
                                        <span class="caret"></span>
                                    </a>
                                    <ul class="dropdown-menu gg-dropd-menu">
                                        <li>
                                            <a class=glygenmenu alt="glycan" href="$moduleDomain/glycan_search.html">Glycan</a>
                                        </li>
                                        <li>
                                            <a class=glygenmenu alt="protein" href="$moduleDomain/protein_search.html">Protein</a>
                                        </li>
                                        <li>
                                            <a class=glygenmenu alt="proteoform" href="$moduleDomain/glycoprotein_search.html">Glycoprotein</a>
                                        </li>
                                    </ul>
                                </li>
                                <li>
                                    <a class=glygenmenu alt="quick search" href="$moduleDomain/index.html#tryMe">TRY&nbsp;ME</a>
                                </li>
                                <li>
                                    <a class=glygenmenu alt="quick search" href="$moduleDomain/quick_search.html">QUICK&nbsp;SEARCH</a>
                                </li>
                                <li>
                                    <a class=glygenmenu alt="about" href="$moduleDomain/about.html">ABOUT</a>
                                </li>
                                
                                <li class="dropdown">
                                    <a class="dropdown-toggle" data-toggle="dropdown" alt="explore" href="$moduleDomain/index.html#">MORE 
                                        <span class="caret"></span>
                                    </a>
                                    <ul class="dropdown-menu gg-dropd-menu">
                                        <li>
                                            <a class=glygenmenu alt="resources" href="$moduleDomain/resources.html">Resources</a>
                                        </li>
                                        <li>
                                            <a class=glygenmenu alt="survey" href="$moduleDomain/survey.html" class="selected">Survey </a>
                                        </li>
                                        <li>
                                            <a class=glygenmenu  alt="contact us" href="$moduleDomain/contact.html">Contact Us</a>
                                        </li>
                                    </ul>
                                </li>
                            </ul>
                            <ul class="nav navbar-nav gg-nav-nav pull-right">
                                <li class="dropdown">
                <!--                    <a href="#" class="dropdown-toggle" data-toggle="dropdown">-->
                                    <a class=glygenmenu href="$moduleDomain/glygen_settings.html">   
                                        <span class="glyphicon glyphicon-user"></span> MY&nbsp;GLYGEN 
                <!--                        <span class="caret"></span>-->
                                    </a>
                <!--
                                    <ul class="dropdown-menu gg-dropd-menu">
                                        <li>
                                            <a class=glygenmenu href="$moduleDomain/glygen_settings.html" class="gg-dropd-backr">
                                                <span class="glyphicon glyphicon-wrench"></span> Privacy Settings
                                            </a>
                                        </li>
                                    </ul>
                -->
                                </li>
                            </ul>
                        </div>
                    </div>
                </nav>
                
                
                </div>
                };
        }

}






#############################
sub getProjectFooter{
        my($moduleDomain, $project) = @_;

        my $libRoot = ($project eq "oncomx" ? qq{$moduleDomain/static/libraries} : qq{$moduleDomain/libraries});        
        if ($project eq "oncomx"){
                return qq{
                <div class="container-fluid text-center">
                <div class="row footer" style="background:#000;">
                <div class="col-sm-12 col-sm-3 copyrights">
                </div>
                <div class="col-sm-12 col-sm-4 policy">
                <ul class="bottom-row-right">
                        <li><a href="$moduleDomain/license.html">License</a></li>
                        <li><a href="$moduleDomain/privacy_policy.html">Privacy&nbsp;Policy</a></li>
                        <li><a href="$moduleDomain/disclaimer.html">Disclaimer</a></li>
                        <li><a href="$moduleDomain/contact.html">Contact&nbsp;Us</a></li>
                </ul>
                </div>
                <div class="col-sm-12 col-sm-2 logos">
                </div>
                </div>
                </div>
                };
        }
        else{
                return qq{

                <div>
                
                    <div class="container-fluid text-center"> 
                    <div class="row footer">
                        <div class="col-sm-12 col-sm-3 copyrights">     
                            &copy; 2018 - GlyGen All rights reserved by <a href="http://www.uga.edu/" target="_blank">UGA</a>&nbsp;and&nbsp;<a href="https://www.gwu.edu/" target="_blank">GWU</a> 
                        </div>  
                
                        <div class="col-sm-12 col-sm-4 policy"> 
                            <ul class="bottom-row-right">
                                <li><a href="$moduleDomain/license.html">License</a></li>
                                <li><a href="$moduleDomain/privacy_policy.html">Privacy&nbsp;Policy</a></li>
                                <li><a href="$moduleDomain/disclaimer.html">Disclaimer</a></li>
                                <li><a href="$moduleDomain/contact.html">Contact&nbsp;Us</a></li>
                            </ul>
                        </div>
                
                        <div class="col-sm-12 col-sm-3 grant">Funded by
                            <a href=" https://commonfund.nih.gov/" target="_blank">NIH Common Funds</a> 
                            <br>
                            <span>Grant # 
                            <a href="http://grantome.com/grant/NIH/U01-GM125267-01" target="_blank">1U01GM125267&nbsp;-&nbsp;01</a>
                            </span>
                        </div> 
                
                        <div class="col-sm-12 col-sm-2 logos">
                            <a href="https://www.ccrc.uga.edu/" target="_blank">
                                <img id="logo-uga" src="/content/logo-uga.png">
                            </a> 
                            <a href="https://smhs.gwu.edu/" target="_blank">
                                <img id="logo-gwu" src="/content/logo-gwu.png">
                            </a>  
                        </div>
                    </div>
                </div>
                
                </div>

                <script type="text/javascript" src="$moduleDomain/libraries/jquery/jquery-3.1.js"></script>
                <script type="text/javascript" src="$moduleDomain/libraries/bootstrap/js/bootstrap.min.js"></script>
                <script src="$moduleDomain/js/navbar.js"></script>
                <script src="$moduleDomain/js/ws_url.js"></script>    
                <!-- script for getting Web Service URLs -->
                <script src="$moduleDomain/js/activity_tracker.js"></script>  
                <!-- script for logging activity -->
                <script src="$moduleDomain/js/utility.js"></script>


                };
        }
}















