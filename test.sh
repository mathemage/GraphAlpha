sudo dpkg --configure -a
find /etc/apt -type f -name '*.list' -exec grep -H '^deb ' '{}' \;
