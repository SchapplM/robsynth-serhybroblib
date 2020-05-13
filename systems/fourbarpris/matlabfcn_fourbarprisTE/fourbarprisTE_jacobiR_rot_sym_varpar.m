% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% JR_rot [9x1]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = fourbarprisTE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_jacobiR_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbarprisTE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_jacobiR_rot_sym_varpar: pkin has to be [3x1] (double)');
JR_rot=NaN(9,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:49
	% EndTime: 2020-05-07 09:01:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:49
	% EndTime: 2020-05-07 09:01:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (152->18), mult. (62->21), div. (16->4), fcn. (8->2), ass. (0->16)
	t50 = qJ(1) + pkin(3);
	t54 = -pkin(2) + t50;
	t46 = pkin(1) + t54;
	t55 = -pkin(2) - t50;
	t47 = pkin(1) + t55;
	t56 = t46 * t47;
	t45 = pkin(1) - t54;
	t53 = t45 * t56;
	t44 = pkin(1) - t55;
	t52 = sqrt(-t44 * t53);
	t51 = 0.1e1 / pkin(1);
	t49 = 1 / (t50 ^ 2);
	t48 = 0.1e1 / t50;
	t43 = (-t50 * t48 - ((-pkin(1) ^ 2 + pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3)) * t49) / 0.2e1) * t51;
	t42 = (t49 * t52 / 0.2e1 - t48 / t52 * (-t53 + (t56 + (t46 - t47) * t45) * t44) / 0.4e1) * t51;
	t1 = [t43; t42; 0; -t42; t43; 0; 0; 0; 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:49
	% EndTime: 2020-05-07 09:01:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (152->18), mult. (62->21), div. (16->4), fcn. (8->2), ass. (0->16)
	t32 = qJ(1) + pkin(3);
	t36 = -pkin(2) + t32;
	t28 = pkin(1) + t36;
	t37 = -pkin(2) - t32;
	t29 = pkin(1) + t37;
	t38 = t28 * t29;
	t27 = pkin(1) - t36;
	t35 = t27 * t38;
	t26 = pkin(1) - t37;
	t34 = sqrt(-t26 * t35);
	t33 = 0.1e1 / pkin(1);
	t31 = 1 / (t32 ^ 2);
	t30 = 0.1e1 / t32;
	t25 = (-t32 * t30 - ((-pkin(1) ^ 2 + pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3)) * t31) / 0.2e1) * t33;
	t24 = (t31 * t34 / 0.2e1 - t30 / t34 * (-t35 + (t38 + (t28 - t29) * t27) * t26) / 0.4e1) * t33;
	t1 = [t24; -t25; 0; 0; 0; 0; t25; t24; 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:49
	% EndTime: 2020-05-07 09:01:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (106->11), mult. (34->10), div. (8->2), fcn. (4->2), ass. (0->13)
	t43 = -pkin(3) - qJ(1);
	t41 = -pkin(2) - t43;
	t35 = pkin(1) + t41;
	t42 = -pkin(2) + t43;
	t36 = pkin(1) + t42;
	t45 = t35 * t36;
	t44 = 0.1e1 / pkin(2) / pkin(1);
	t34 = pkin(1) - t41;
	t40 = t34 * t45;
	t33 = pkin(1) - t42;
	t39 = (-t33 * t40) ^ (-0.1e1 / 0.2e1) * (-t40 + (t45 + (t35 - t36) * t34) * t33) * t44;
	t32 = t43 * t44;
	t1 = [t32; -t39 / 0.4e1; 0; t39 / 0.4e1; t32; 0; 0; 0; 0;];
	JR_rot = t1;
end