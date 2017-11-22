// Open all the image files of a given extension in the selected folder

macro "Open Folder" {

  // Get the folder to open
  input = getDirectory("Input directory");

  // Get the contained files
  list = getFileList(input);

  // Get the file extension to open
  Dialog.create("File type");
  Dialog.addString("File suffix: ", ".tif", 5);
  Dialog.show();
  suffix = Dialog.getString();

  // Open all files ...
  for (i = 0; i < list.length; i++) {

    // ... that have the corresponding extension
    if(endsWith(list[i], suffix)) {
      file = list[i];

      // Print the name
      print("Opening: " + input + file);
      open(input + file);
    }
  }
}
