bl_info = {
	"name": "(FLIM) Fluorescent Life-time Imaging Microscopy",
	"description": "Import/Export *.flim files obtained from processing the measurements",
	"author": "Denys I. Bondar",
	"version": (1, 0),
	"blender": (2, 70, 0),
	"location": "File > Import > Import FLIM data and File > Export > Export FLIM data ",
	"warning": "", # used for warning icon and text in addons panel
	#"wiki_url": "http://wiki.blender.org/index.php/Extensions:2.5/Py/Scripts/My_Script",
	"category": "Import-Export"
}

import bpy, pickle, operator, hashlib

def get_file_hash_sum (filepath) :
	"""
	Generate hashume from file specified by its name.
	"""
	m = hashlib.sha1()
	with open(filepath, "rb") as in_file :
		for chunk in in_file : m.update(chunk)
	return m.hexdigest()
	
def import_flim_data (filepath) :
	"""
	Importing FLIM files (*.flim) created by the method <save_post_processed> in <viewer.py>
	"""
	# loading data
	with open(filepath, "rb") as in_file : fitted_data = pickle.load(in_file) 
		
	# Creating a group
	FLIM_group = bpy.data.groups.new("FLIM_group")
	# Saving the file name and the file hash sum 
	FLIM_group["film_import_filename"] = filepath
	FLIM_group["film_import_file_hash_sum"] = get_file_hash_sum(filepath)
	
	print ("\nFLIM file successfully loaded\n")
	
	# Extract the conversion factors from units to microns
	unit_to_microns_ax1 = fitted_data["unit_to_microns_ax1"]
	unit_to_microns_ax2 = fitted_data["unit_to_microns_ax2"]
	unit_to_microns_ax3 = fitted_data["unit_to_microns_ax3"]
	
	cm_ax1, cm_ax2, cm_ax3 = fitted_data["centre_of_mass"]
	
	# pixel size in microns
	pixel_size = fitted_data["pixel_size"]
	pixel_size = ( 0.52*pixel_size[0]*unit_to_microns_ax1, 0.52*pixel_size[1]*unit_to_microns_ax2, 0.52*pixel_size[2]*unit_to_microns_ax3 )
	
	print ("\nImporting materials...\n")
	
	# Adding materials
	for key, material in fitted_data["materials"].items() :
		bpy.data.materials.new (name=key)
		M = bpy.data.materials[key]
		M.use_transparency = True
		for property_name, property_value in material.items() : 
			setattr(M, property_name, property_value)
		 
	print ("\nImporting pixels...\n")
	
	# Functions for fast access 
	get_material = operator.itemgetter("material")
	
	# To improve performance in Blender, we create a pixel then we copy it 
	bpy.ops.mesh.primitive_cube_add ()
	bpy.ops.transform.resize (value=pixel_size)
	original_pixel = bpy.context.object
	
	# Looping over points
	for key, value in fitted_data["points"].items() :
		# Get location in microns
		location = ( unit_to_microns_ax1*(key[0]-cm_ax1), unit_to_microns_ax2*(key[1]-cm_ax2), unit_to_microns_ax3*(key[2]-cm_ax3) )
		
		# Copying 
		current_pixel = original_pixel.copy()
		current_pixel.data = original_pixel.data.copy()
		current_pixel.location = location
		
		# Changing the attributes of the current pixel
		current_pixel.data.materials.append (bpy.data.materials[get_material(value)])
		current_pixel.name = "pixel_%d_%d_%d" % key
		
		# Saving key and value
		current_pixel.data["FLIM_key"] = key
		current_pixel.data["FLIM_value"] = value
		
		# Adding pixel to scene
		bpy.context.scene.objects.link(current_pixel)
		# Adding object to grope
		FLIM_group.objects.link (current_pixel)
	
	bpy.context.scene.update()
	print ("\nFLIM file imported!!!!\n")
	
def export_flim_data (filepath) :
	"""
	Export FLIM file 
	"""
	# Re-loading data from the imported file
	film_import_filename = bpy.data.groups["FLIM_group"]["film_import_filename"]
	with open(film_import_filename, "rb") as in_file :
		fitted_data = pickle.load(in_file) 
	 
	# Verifying hash-sum to make sure its the same file
	if bpy.data.groups["FLIM_group"]["film_import_file_hash_sum"] != get_file_hash_sum(film_import_filename) :
		raise IOError ("Export failed. Because imported file <%s> has been modified! " % film_import_filename)

	print("\nImported file reloaded.\n")
	
	# Deleting the info about points
	del fitted_data["points"]
	
	# Functions for fast access 
	get_FLIM_key 	= operator.itemgetter("FLIM_key")
	get_FLIM_value 	= operator.itemgetter("FLIM_value")
	
	# Collecting information about points
	points = {}
	for obj in bpy.data.objects :
		try : 
			points[ tuple(get_FLIM_key(obj.data)) ] = get_FLIM_value(obj.data).to_dict()
		except KeyError : pass
	
	fitted_data["points"] = points
	
	# Saving the file
	with open(filepath, "wb") as out_file : 
		# using backward compatible protocol
		pickle.dump(fitted_data, out_file, protocol=2) 
	
	# Updating the hash-sum since the file has been edited 
	bpy.data.groups["FLIM_group"]["film_import_file_hash_sum"] = get_file_hash_sum(filepath)

	print("\nExport successfully completed\n")


class CExport_FLIM_Data(bpy.types.Operator):
	"""
	Class for exporting FLIM measurements
	"""
	bl_idname = "export.flim_data"
	bl_label = "Export FLIM Data"
	filepath = bpy.props.StringProperty(subtype="FILE_PATH")
	filter_glob = bpy.props.StringProperty(default="*.flim", options={'HIDDEN'})
	
	def execute(self, context):
		export_flim_data(self.filepath)
		return {'FINISHED'}

	def invoke(self, context, event):
		context.window_manager.fileselect_add(self)
		return {'RUNNING_MODAL'}

class CImport_FLIM_Data(bpy.types.Operator):
	"""
	Class for importing FLIM measurements
	"""
	bl_idname = "import.flim_data"
	bl_label = "Import FLIM Data"
	filepath = bpy.props.StringProperty(subtype="FILE_PATH")
	filter_glob = bpy.props.StringProperty(default="*.flim", options={'HIDDEN'})
		
	def execute(self, context):
		import_flim_data(self.filepath)
		return {'FINISHED'}

	def invoke(self, context, event):
		context.window_manager.fileselect_add(self)
		return {'RUNNING_MODAL'}
		

# Only needed if you want to add into a dynamic menu
def export_menu_func(self, context):
	self.layout.operator_context = 'INVOKE_DEFAULT'
	self.layout.operator(CExport_FLIM_Data.bl_idname, text="Export FLIM data (*.flim)...")

def import_menu_func(self, context):	
	self.layout.operator_context = 'INVOKE_DEFAULT'
	self.layout.operator(CImport_FLIM_Data.bl_idname, text="Import FLIM data (*.flim)...")

# register the class
def register():
	#bpy.utils.register_module(__name__)
	bpy.utils.register_class (CExport_FLIM_Data)
	bpy.utils.register_class (CImport_FLIM_Data)
	
	# add a menu item into the file selector
	bpy.types.INFO_MT_file_export.append(export_menu_func)
	bpy.types.INFO_MT_file_import.append(import_menu_func)


def unregister():
	#bpy.utils.unregister_module(__name__)
	bpy.utils.unregister_class (CExport_FLIM_Data)
	bpy.utils.unregister_class (CImport_FLIM_Data)
	
	# remove the menu item
	bpy.types.INFO_MT_file_export.remove(export_menu_func)
	bpy.types.INFO_MT_file_import.remove(import_menu_func)

if __name__ == "__main__":
	register()
