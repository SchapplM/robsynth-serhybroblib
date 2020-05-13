% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh1m2DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh1m2DE2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_jacobia_rot_sym_varpar: pkin has to be [22x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:33
	% DurationCPUTime: 1.21s
	% Computational Cost: add. (10672->25), mult. (18186->50), div. (211->9), fcn. (27825->20), ass. (0->43)
	t146 = pkin(22) + pkin(21);
	t144 = sin(t146);
	t145 = cos(t146);
	t147 = cos(pkin(20));
	t152 = sin(pkin(18));
	t180 = sin(pkin(20));
	t181 = cos(pkin(18));
	t141 = t152 * t147 - t181 * t180;
	t142 = t181 * t147 + t152 * t180;
	t149 = sin(qJ(3));
	t154 = cos(qJ(3));
	t136 = t141 * t154 - t149 * t142;
	t150 = sin(qJ(2));
	t155 = cos(qJ(2));
	t183 = t149 * t141 + t142 * t154;
	t182 = -t136 * t155 + t150 * t183;
	t188 = t150 * t136 + t183 * t155;
	t124 = qJ(2) + qJ(3) + atan2(t188 * t144 + t182 * t145, t182 * t144 - t145 * t188);
	t122 = sin(t124);
	t119 = t122 ^ 2;
	t123 = cos(t124);
	t151 = sin(qJ(1));
	t114 = 0.1e1 / (t151 ^ 2 * t119 / t123 ^ 2 + 0.1e1);
	t195 = t114 / t123;
	t166 = t151 * t122;
	t113 = atan2(-t166, -t123);
	t111 = sin(t113);
	t112 = cos(t113);
	t107 = -t111 * t166 - t112 * t123;
	t156 = cos(qJ(1));
	t177 = 0.1e1 / t107 ^ 2 * t156 ^ 2;
	t148 = sin(qJ(4));
	t165 = t151 * t148;
	t153 = cos(qJ(4));
	t164 = t151 * t153;
	t163 = t156 * t148;
	t162 = t156 * t153;
	t118 = t123 * t162 + t165;
	t116 = 0.1e1 / t118 ^ 2;
	t117 = t123 * t163 - t164;
	t161 = t117 ^ 2 * t116 + 0.1e1;
	t110 = 0.1e1 / t161;
	t1 = [t122 * t156 * t195, 0, 0, 0; (-0.1e1 / t107 * t166 - (-t112 * t119 * t151 * t195 + (t114 - 0.1e1) * t122 * t111) * t122 * t177) / (t119 * t177 + 0.1e1), 0, 0, 0; ((-t123 * t165 - t162) / t118 - (-t123 * t164 + t163) * t116 * t117) * t110, 0, 0, t161 * t110;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end