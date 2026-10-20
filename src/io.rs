// Copyright 2026 Henrik Mengkrogen
// Developed for the Aarhus University iGEM team.
//
// Licensed under the Apache License, Version 2.0.
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use std::fs::{self};
use std::io::{self, Write};

use std::process::{Command, Stdio};



pub fn read_input_file(path: &str) -> Result<(String, String), std::io::Error> {
    let content = fs::read_to_string(path)?;

    let mut seq = String::new();
    let mut target = String::new();

    for line in content.lines() {
        let line = line.trim();

        if line.starts_with("Sequence") {
            if let Some(idx) = line.find(':') {
                seq = line[idx + 1..].trim().to_string();
            }
        } else if line.starts_with("Structure") {
            if let Some(idx) = line.find(':') {
                target = line[idx + 1..].trim().to_string();
            }
        }
    }

    Ok((seq, target))
}

pub fn ask_positive_usize(prompt: &str, default: usize) -> usize {
    loop {
        print!("{prompt} [{default}]: ");
        io::stdout()
            .flush()
            .expect("Failed to flush terminal output");

        let mut input = String::new();
        io::stdin()
            .read_line(&mut input)
            .expect("Failed to read terminal input");

        let input = input.trim();

       
        if input.is_empty() {
            return default;
        }

        match input.parse::<usize>() {
            Ok(value) if value > 0 => return value,
            _ => {
                println!("Please enter a positive whole number.");
            }
        }
    }
}

pub fn ask_positive_i64(prompt: &str, default: i64) -> i64 {
    loop {
        print!("{prompt} [{default}]: ");
        io::stdout()
            .flush()
            .expect("Failed to flush terminal output");

        let mut input = String::new();
        io::stdin()
            .read_line(&mut input)
            .expect("Failed to read terminal input");

        let input = input.trim();

        // Pressing Enter accepts the default.
        if input.is_empty() {
            return default;
        }

        match input.parse::<i64>() {
            Ok(value) if value > 0 => return value,
            _ => {
                println!("Please enter a positive whole number.");
            }
        }
    }
}

pub fn ask_to_view_results() -> io::Result<bool> {
    print!("\nView all final results in a scrollable window? (y/n): ");
    io::stdout().flush()?;

    let mut answer = String::new();
    io::stdin().read_line(&mut answer)?;

    Ok(answer.trim().eq_ignore_ascii_case("y"))
}

pub fn show_in_pager(output: &str) -> io::Result<()> {
    let mut pager = Command::new("less")
        .arg("-R")
        .stdin(Stdio::piped())
        .spawn()?;

    if let Some(mut pager_input) = pager.stdin.take() {
        pager_input.write_all(output.as_bytes())?;
    }

    pager.wait()?;
    Ok(())
}

pub fn gc_content(sequence: &str) -> Option<f64> {
    let mut gc_count = 0usize;
    let mut base_count = 0usize;

    for base in sequence.bytes() {
        match base.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc_count += 1;
                base_count += 1;
            }
            b'A' | b'U' => {
                base_count += 1;
            }
            _ => {}
        }
    }

    if base_count == 0 {
        return None;
    }

    let percentage = (gc_count as f64 / base_count as f64) * 100.0;
    Some((percentage * 100.0).round() / 100.0)
}

