% Jacobian of implicit kinematic constraints of
% fivebar1IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = fivebar1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:08
% EndTime: 2020-04-27 06:19:08
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->6), mult. (8->8), div. (0->0), fcn. (8->8), ass. (0->3)
t4 = qJ(1) + qJ(2);
t3 = qJ(3) + qJ(4);
t1 = [pkin(5) * sin(t4) - sin(qJ(1)) * pkin(2), pkin(4) * sin(t3) + sin(qJ(3)) * pkin(3); -pkin(5) * cos(t4) + cos(qJ(1)) * pkin(2), -pkin(4) * cos(t3) - cos(qJ(3)) * pkin(3); 1, -1;];
Phi_a = t1;
