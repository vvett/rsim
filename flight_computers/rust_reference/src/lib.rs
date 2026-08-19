use std::ffi::{c_char, c_void};

#[repr(C)]
pub struct Descriptor {
    struct_size: u32,
    abi_version: u32,
    plugin_name: *const c_char,
    sensor_count: usize,
    sensor_names: *const *const c_char,
    actuator_count: usize,
    actuator_names: *const *const c_char,
    configuration_count: usize,
    configuration_names: *const *const c_char,
    telemetry_count: usize,
    telemetry_names: *const *const c_char,
}

#[repr(C)]
pub struct Configuration {
    struct_size: u32,
    name: *const c_char,
    value: f64,
}

unsafe impl Sync for Descriptor {}

static PLUGIN_NAME: &[u8] = b"rust_reference\0";
static DESCRIPTOR: Descriptor = Descriptor {
    struct_size: std::mem::size_of::<Descriptor>() as u32,
    abi_version: 2,
    plugin_name: PLUGIN_NAME.as_ptr() as *const c_char,
    sensor_count: 0,
    sensor_names: std::ptr::null(),
    actuator_count: 0,
    actuator_names: std::ptr::null(),
    configuration_count: 0,
    configuration_names: std::ptr::null(),
    telemetry_count: 0,
    telemetry_names: std::ptr::null(),
};

#[no_mangle]
pub extern "C" fn rsim_flight_computer_get_descriptor() -> *const Descriptor {
    &DESCRIPTOR
}

#[no_mangle]
pub extern "C" fn rsim_flight_computer_create(
    _configuration: *const Configuration,
    configuration_count: usize,
    _error: *mut c_char,
    _error_size: usize,
) -> *mut c_void {
    if configuration_count != 0 {
        return std::ptr::null_mut();
    }
    Box::into_raw(Box::new(0_u8)) as *mut c_void
}

#[no_mangle]
pub unsafe extern "C" fn rsim_flight_computer_step(
    _state: *mut c_void,
    _time: f64,
    _dt: f64,
    _inputs: *const f64,
    input_count: usize,
    _outputs: *mut f64,
    output_count: usize,
    _telemetry: *mut f64,
    telemetry_count: usize,
    _error: *mut c_char,
    _error_size: usize,
) -> i32 {
    if input_count == 0 && output_count == 0 && telemetry_count == 0 {
        0
    } else {
        1
    }
}

#[no_mangle]
pub unsafe extern "C" fn rsim_flight_computer_destroy(state: *mut c_void) {
    if !state.is_null() {
        drop(Box::from_raw(state as *mut u8));
    }
}
