/*
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2020.
 *
 *  @file  SpatialRestraint-web.cc
 *  @brief Web app to show the growth of a multicell under different conditions
 *  @note Status: BETA
 */


// TODO: Make colorblind friendly
// TODO: Add a description of what's going on
// TODO: Add color key
// TODO: Look into making updates faster!
  // If unrestrained cells are present, we may be redrawing the same cell many times
    // But is this slower than redrawing all cells every time? 

// Standard

// Empirical
#include "emp/web/web.hpp"
#include "emp/web/Animate.hpp"
#include "emp/web/color_map.hpp"
#include "emp/web/Element.hpp"

// Local
#include "../SpatialRestraint.h"

class RogueCellAnimationController: public emp::web::Animate{
private:
  emp::web::Document doc;                   // Div that all our elements fall in
  const int canvas_width = 512;             // Width in pixels
  const int canvas_height = 512;            // Height in pixels
  int mc_size = 128;                        // Number of cells on one side of the multicell
  int tile_width = canvas_width / mc_size;  // Width of each cell in pixels
  int tile_height = canvas_height / mc_size;// Height of each cell in pixels
  int steps_per_draw = 1024;                // How many cell births per canvas update?
  size_t starting_ones = 55;                // The number of ones in the starting cell's genome
  double mut_prob = 0.2;                    // Per genome mutation rate (only ever one bit flip)
  bool one_check = false;                   // If true, restrained cells only check one random 
                                            //    neighbor before losing their resources. If false, 
                                            //    restrained cells select randomly from empty 
                                            //    neighbors
  bool run_to_end = false;                  // Do we render multiple updates or only the very end?
  Multicell multicell;                      // The multicell to be visualized
  emp::Random random;                       // Random number generator
  emp::vector<size_t> cells_to_draw;        // List of cell ids to be rendered
  size_t cells_to_draw_count;               // How many cells have been placed / need rendered?

protected:
  // Takes the index of a cell in the multicell, and draws it to the canvas
  void DrawCell(size_t cell_id){
    auto canvas = doc.Canvas("canvas");
    size_t num_ones = multicell.cells[cell_id].num_ones;
    std::string color_fill = "#ffffff";
    if(num_ones < 50){ // Unrestrained
        color_fill = emp::ColorRGB(
          255 - (49 - num_ones) * 5,
          0,
          0 + (49 - num_ones) * 2);
    }
    else{ // Restrained
        color_fill = emp::ColorRGB(
          255 - (num_ones - 50) * 5,
          255 - (num_ones - 50) * 5,
          255 - (num_ones - 50) * 5);
    }
    canvas.Rect(
      (cell_id % mc_size) * tile_width, 
      (cell_id / mc_size) * tile_height,
      tile_width, 
      tile_height,
      color_fill,
      "#00000000"); // No edge lines!
  }
  // Called for every update of the animation. 
  // If run_to_end is true, simulate the multicell until full, then draw
  // Else, simulate steps_per_draw attempted cell reproductions and draw
  void DoFrame() override {
    if(run_to_end){ 
      while(multicell.num_cells < multicell.cells.size()){
          multicell.DoStep();
      }
      for(size_t idx = 0; idx < mc_size * mc_size; ++idx){
        DrawCell(idx);
      }
      Stop();
      doc.Button("anim_toggle_btn").SetLabel("Start");
      doc.Button("anim_toggle_btn").SetAttr("class", "btn btn-success");
    }
    else{
      cells_to_draw_count = 0;
      if(multicell.num_cells < multicell.cells.size()){
        int actual_steps = steps_per_draw;
        if(steps_per_draw == -1){
          actual_steps = multicell.num_cells;
        }
        for(size_t step_id = 0; step_id < actual_steps; ++step_id){
          multicell.DoStep();
          if(multicell.cell_placed_last_step){
            cells_to_draw[cells_to_draw_count++] = multicell.last_placed_cell_id;
          }
        }
        for(size_t idx = 0; idx < cells_to_draw_count; ++idx){
          DrawCell(cells_to_draw[idx]);
        }
      }
      else{ // We've finished! Stop the animation and update the buttons
        Stop();
        doc.Button("anim_toggle_btn").SetLabel("Start");
        doc.Button("anim_toggle_btn").SetAttr("class", "btn btn-success");
      }
    }
  }
  // Reset the simulation. Stop animation, restart the multicell from a single cell.
  void Reset(){
    Stop();
    doc.Button("anim_toggle_btn").SetLabel("Start");
    doc.Button("anim_toggle_btn").SetAttr("class", "btn btn-success");
    multicell.SetupConfig();
    auto canvas = doc.Canvas("canvas");
    canvas.Clear("#000000");
    multicell.InjectCell(mc_size / 2 +  mc_size * (mc_size / 2), starting_ones);
    cells_to_draw.resize(0);
    DrawCell(mc_size / 2 + mc_size * (mc_size / 2));
  }
  // Pull configuration options from the web page and apply them to the multicell
  void ApplyConfig(){
    mc_size = std::stoi(doc.Input("mc_size_input").GetCurrValue());
    tile_width = canvas_width / mc_size;
    tile_height = canvas_height / mc_size;
    starting_ones = std::stoi(doc.Input("ones_input").GetCurrValue());
    mut_prob = std::stod(doc.Input("mut_prob_input").GetCurrValue());
    std::string one_check_str = doc.Input("one_check_input").GetCurrValue();
    one_check = (one_check_str == "1" || one_check_str == "true");
    multicell.mut_prob = mut_prob;
    multicell.cells_side = mc_size;
    multicell.restrain = 50;
    multicell.genome_size = 100;
    multicell.one_check = one_check;
    multicell.SetupConfig();
  }
  // Create all the infrastructure to have an input with Bootstrap's input-group and add-on text
  // Returns div of input_group
  emp::web::Div AddInputGroup(const std::string & id_prefix, const std::string & addon_text, 
      const std::string & input_type, const std::function<void(const std::string &)> &  callback){
    auto input_group = emp::web::Div(id_prefix + "_input_group").SetAttr("class", "input-group");
    auto addon = emp::web::Div(id_prefix + "_addon_text").SetAttr("class", "input-group-addon");
    addon << addon_text;
    input_group << addon;
    auto input = emp::web::Input(callback, input_type, "", id_prefix + "_input");
    input.SetCSS("width", "99%");
    input_group << input;
    return input_group;
  }
  // Add widgets to allow the user to change configuration options
  void AddConfigWidgets(){
    auto config_col = doc.Div("col_config"); 
    auto config_panel = emp::web::Div("config_panel").SetAttr("class", "panel panel-default");
    config_col << config_panel;
    // Header text
    auto config_panel_header =emp::web::Div("config_panel_header").SetAttr("class", "panel-heading");
    config_panel << config_panel_header;
    config_panel_header << "<center><h3>Configuration</h3></center>";
    // Meat of the configuration
    auto config_panel_body = emp::web::Div("config_panel_body").SetAttr("class", "panel-body");
    config_panel << config_panel_body;
    auto form = emp::web::Element("form", "config_form");
    config_panel_body << form;
    // MC size
    auto mc_size_input_group = AddInputGroup("mc_size", "Cells per side", "number", 
      [](const std::string & s){;});
    form << mc_size_input_group; 
    doc.Input("mc_size_input").Value(mc_size);
    // Starting ones
    auto ones_input_group = AddInputGroup("ones", "Ones in genome at start", "number", 
      [](const std::string & s){;});
    form << ones_input_group; 
    doc.Input("ones_input").Value(starting_ones);
    // Mutation rate
    auto mut_prob_input_group = AddInputGroup("mut_prob", "Mutation probability", "number", 
      [](const std::string & s){;});
    form << mut_prob_input_group; 
    doc.Input("mut_prob_input").Value(mut_prob);
    // One check (i.e., do restrained cells only check on neighbor at random?)
    auto one_check_input_group = AddInputGroup("one_check",
        "Restrained cells check only a single neighbor?", "checkbox", [](const std::string & s){;});
    form << one_check_input_group; 
    doc.Input("one_check_input").Value(one_check);
    config_panel_body << "<br/>";
    // Apply button
    auto center = emp::web::Element("center", "");
    config_panel_body << center;
    center << emp::web::Button(
      [this](){ Stop(); ApplyConfig(); Reset(); }, 
      "Restart & Apply", "config_apply_btn").SetAttr("class", "btn btn-primary");
  }
  // Add the canvas and the widgets that directly control it. 
  void AddCanvas(){
    // DOM Strucutre
    //  - Canvas bootstrap column
    //    - Canvas Panel
    //      - Canvas Panel Header
    //      - Canvas Panel Body
    //        - Canvas (the actualy draw area
    //        - Canvas controls (toggle, step, reset, step selection)
    auto canvas_col = doc.Div("col_canvas");
    auto canvas_panel = emp::web::Div("canvas_panel").SetAttr("class", "panel panel-default");
    canvas_col << canvas_panel;
    auto canvas_panel_header =emp::web::Div("canvas_panel_header").SetAttr("class", "panel-heading");
    canvas_panel << canvas_panel_header;
    canvas_panel_header << "<center><h3>Visualization</h3></center>";
    auto canvas_panel_body = emp::web::Div("canvas_panel_body").SetAttr("class", "panel-body");
    canvas_panel << canvas_panel_body;
    auto canvas = emp::web::Canvas(canvas_width, canvas_height, "canvas");
    targets.push_back(canvas);
    canvas.SetCSS("display", "block");
    canvas.SetCSS("margin", "0 auto");
    canvas_panel_body << canvas;
    canvas.Clear("#000000");
    canvas_panel_body << "<br/>";
    auto center = emp::web::Element("center", "");
    canvas_panel_body << center;
    auto canvas_toggle_button = emp::web::Button([this](){
        ToggleActive();
        doc.Button("anim_toggle_btn").SetLabel(active ? "Stop" : "Start");
        doc.Button("anim_toggle_btn")
          .SetAttr("class", (active ? "btn btn-danger" : "btn btn-success"));
        }, "Start", "anim_toggle_btn").SetAttr("class", "btn btn-success");
    center << canvas_toggle_button;
    //center << GetToggleButton("anim_toggle_btn", "Start");
    center << "&nbsp";
    center << GetStepButton("anim_step_btn", "Advance Step")
      .SetAttr("class", "btn btn-default");
    center << "&nbsp";
    center << emp::web::Button([this](){ Reset(); }, "Reset Multicell", "anim_reset_btn")
      .SetAttr("class", "btn btn-warning");
    canvas_panel_body <<  "<br/>";
    // Allow user to change the number of attempted cell births per render update
    auto form = emp::web::Element("form", "canvas_form");
    canvas_panel_body << form;
    auto steps_input_group = AddInputGroup("steps", "Attempted cell births per draw (-1 for auto)",
      "number", [this](const std::string & s){ 
        steps_per_draw = std::stoi(s);
    }); 
    form << steps_input_group; 
    doc.Input("steps_input").Value(steps_per_draw);
    // Should we just simulate all the way to the end as one update?
    auto run_to_end_input_group = AddInputGroup("run_to_end",
      "Jump to end?", "checkbox", [this](const std::string & s){
          run_to_end = (s == "1" || s == "true");
      });
    form << run_to_end_input_group; 
    doc.Input("run_to_end_input").Value(run_to_end);
  }
public:
  RogueCellAnimationController():
    doc("emp_base"),
    multicell(random)
  {
    cells_to_draw.resize(512 * 512, 0); // Largest possible size
    // Setup the basics of bootstrap (overall container, jumbotron, etc.)
    auto container_div = doc.AddDiv("container_main").SetAttr("class", "container");
    auto header = emp::web::Div("header").SetAttr("class", "jumbotron");
    header << "<center><h1>Rogue Cell Multicell Visualization</h1></center>";
    container_div << header;
    auto bootstrap_row = emp::web::Div("row_main").SetAttr("class", "row");
    container_div << bootstrap_row;
    auto canvas_col = emp::web::Div("col_canvas").SetAttr("class", "col-md-6");
    bootstrap_row << canvas_col;
    bootstrap_row << emp::web::Div("col_config").SetAttr("class", "col-md-6");
    // Add all the interactive bits
    AddCanvas(); 
    AddConfigWidgets();
    // Reset all configuration options
    ApplyConfig();
    multicell.SetupConfig();
    // Insert and draw the first (center) cell!
    multicell.InjectCell(mc_size / 2 +  mc_size * (mc_size / 2), starting_ones);
    DrawCell(mc_size / 2 + mc_size * (mc_size / 2));
  }

};

RogueCellAnimationController anim_controller;
int main(){
}

void empCppCallback(){
}
