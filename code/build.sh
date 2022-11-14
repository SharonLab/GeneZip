#!/usr/bin/env bash

if ! command -v cargo &> /dev/null; then
    script_path=$( mktemp )
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs > "${script_path}"
    sh "${script_path}" -qy
    rm "${script_path}"
    source "$HOME/.cargo/env"
    rustup toolchain install stable-x86_64-unknown-linux-gnu
fi

cargo build --release
