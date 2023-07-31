# Return the name of the pair referenced in a file path
def pull_pair_from_string(s):
    strings_to_check = ["piplus_piplus", "piplus_pi0", "piplus_piminus", "pi0_pi0", "piminus_pi0", "piminus_piminus"]

    for string in strings_to_check:
        if string in s:
            return string

    return None

# Given a pion pair string (ex: piplus_piminus) return a string compatible with LaTeX
def latex_pion_pair(pair,use_hashtag=False):
    if use_hashtag:
        return pair.replace("_","").replace("piplus","#pi^{+}").replace("piminus","#pi^{-}").replace("pi0","#pi^{0}")
    else:
        return pair.replace("_","").replace("piplus","\pi^{+}").replace("piminus","\pi^{-}").replace("pi0","\pi^{0}")   
    
# Return the version (aka run period) referenced in the string
def pull_version_from_string(s):
    strings_to_check = ["Fall2018_RGA_inbending", "Fall2018_RGA_outbending","Fall2018Spring2019_RGA_inbending","Spring2019_RGA_inbending"]

    for string in strings_to_check:
        if string in s:
            return string

    return None

# Return a compactified run period
def printable_version(version):
    version = version.replace("Fall2018Spring2019_RGA_inbending","inb. rg-a")

    version = version.replace("Fall2018_RGA_inbending","inb. f18 rg-a")
    version = version.replace("Fall2018_RGA_outbending","outb. f18 rg-a")
    version = version.replace("Spring2019_RGA_inbending","inb. sp19 rg-a")
    version = version.replace("MC_RGA_inbending","inb. MC rg-a")
    version = version.replace("MC_RGA_outbending","outb. MC rg-a")

    version = version.replace("Spring2020_RGB_inbending","inb. sp20 rg-b")
    version = version.replace("Fall2019_RGB_outbending","outb. f19 rg-b")
    version = version.replace("Spring2019_RGB_inbending","inb. sp19 rg-b")
    version = version.replace("MC_RGB_inbending","inb. MC rg-b")
    version = version.replace("MC_RGB_outbending","outb. MC rg-b")
    
    return version