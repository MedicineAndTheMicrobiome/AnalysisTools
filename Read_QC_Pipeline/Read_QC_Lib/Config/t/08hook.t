#
# Testcases for using perl hooks in Ini configuration files
# GCARLS 04/29/2005
#

use strict;
use Test;

BEGIN { $| = 1; plan tests => 13 }
use Config::IniFiles;
my $loaded = 1;
#Test 1
ok($loaded);

my $ini;

# Get files from the 't' directory, portably
chdir('t') if ( -d 't' );

# Test 2
# Create ini object with allowcode Option set to 1
$ini = new Config::IniFiles(-file => 'hook.ini', 
 		            -allowcode => 1);
ok($ini->{allowcode});

#Test 3
# check standard parameter access
ok($ini->val('hooksection', 'testval') eq 'ok');

#Test 4
#check perl hook - accessing environment variables from a hooks sub
$ENV{'HOOK'}='PERLHOOK';
ok($ini->val('hooksection', 'envhook') eq 'PERLHOOK');

#Test 5
# check if a global sub can be called from a hook

sub gethook {
  my($arg)=@_;
  return("HOOK($arg)");
}
ok($ini->val('hooksection', 'hooksub') eq 'HOOK(4711)');

#Test 6
#check if single elements of an array are evaluated
my(@hookary);
@hookary=$ini->val('hooksection', 'hookary');
ok($#hookary==3);

#Test 7
ok($hookary[0] eq 'hook1');

#Test 8
ok($hookary[1] eq 'PERLHOOK');

#Test 9
ok($hookary[2] eq 'HOOK(4711)');

#Test 10
ok($hookary[3] eq 'hook4');

#Test 11
# Rewrite File
$ini->WriteConfig('hook_2.ini');

$ENV{'HOOK'}='PERLHOOK_2';
# check if perl hook still exist in the written file
$ini = new Config::IniFiles(-file => 'hook_2.ini', 
 		            -allowcode => 1);
ok($ini->val('hooksection', 'envhook') eq 'PERLHOOK_2');

#Test 12
# check if perl hook in here document still exists
@hookary=$ini->val('hooksection', 'hookary');
ok($hookary[1] eq 'PERLHOOK_2');


#Test 13
# check if -allowcode => 0 prohibits perl code in ini files
$ini = new Config::IniFiles(-file => 'hook.ini', 
 		            -allowcode => 0);

eval{
  $ini->val('hooksection', 'hooksub');
  ok(0);
} or ok(1);

unlink 'hook_2.ini';
