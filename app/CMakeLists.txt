# SPDX-FileCopyrightInfo: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dumux_add_test(NAME test_erosion SOURCES main.cc)
